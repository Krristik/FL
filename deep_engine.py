import random
import os
import json
import subprocess
import tempfile
import time
import math
from typing import List, Optional
from deap import base, tools
import csv

import requests
from scipy.constants import Boltzmann

from shared.fed_node.node_state import NodeState
from shared.logging_config import logger
from fog_node.fog_resources_paths import FogResourcesPaths
from shared.shared_resources_paths import SharedResourcesPaths
from fog_node.genetics.genetic_engine import Individual, FitnessMulti, select_evaluation_node
from shared.fed_node.fed_node import MessageScope
from shared.utils import metric_weights, reinitialize_and_set_parent


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

DEEP_CSV_LOG_PATH = "/app/cache/deep_params_log.csv"
DEEP_BINARY_PATH = "/usr/local/bin/deepmethod"
DEEP_WRAPPER_PATH = "/app/fog_node/deep_wrapper_logger.py"


class DeepEngine:
    """
    DeepEngine replaces GeneticEngine by delegating hyperparameter optimisation
    to the DEEP (Differential Evolution Entirely Parallel) binary.

    How it works per FL round:
      1. evolve() clears the CSV log and selects an evaluation edge node.
      2. A DEEP config is generated (bounds, algorithm settings).
      3. Node-specific variables (MODEL_FILE_PATH, VALIDATION_NODE_URL, etc.)
         are injected via subprocess env — [default_env] is NOT used because
         this DEEP build ignores that section.
      4. DEEP calls deep_wrapper_logger.py for each individual evaluation.
      5. The wrapper POSTs to the edge node and logs (hyperparams, MSE) to CSV.
      6. After DEEP exits, evolve() reads the CSV and builds a ranked population.

    Config notes (verified against working manual run):
      - type=command  (not type=external — causes GLib crash)
      - Semicolon delimiters on all list fields in [default_model]
      - es_lambda=2   (higher values caused g_async_queue assertion crash)

    Current limitations (intentional):
      - Only MSE (after_training) is used as the fitness objective.
      - resource_load is always 0.0.
      - stagnation_counter is updated for interface compatibility only.
    """

    def __init__(self):
        self.population_size = 5
        self.number_of_generations = 3
        self.stagnation_limit = 2

        self.toolbox = None
        self.best_fitness = float("inf")
        self.stagnation_counter = 0
        self.boltzmann_constant = Boltzmann
        self.additional_factor = 1e23
        self.current_population = None
        self.operating_data_date = []
        self.genetic_evaluation_strategy = None
        self.fog_rabbitmq_host = None
        self.fog_edge_send_exchange = None
        self.edge_to_fog_queue = None

        # kept for interface compatibility; not adapted at runtime
        self.crossover_probability = 0.7
        self.mutation_probability = 0.2
        self.min_crossover_probability = 0.3
        self.max_crossover_probability = 0.9
        self.min_mutation_probability = 0.1
        self.max_mutation_probability = 0.6

        self.fog_model_file_path = os.path.join(
            FogResourcesPaths.MODELS_FOLDER_PATH,
            FogResourcesPaths.FOG_MODEL_FILE_NAME
        )
        self.genetic_population_file_path = os.path.join(
            SharedResourcesPaths.CACHE_FOLDER_PATH,
            FogResourcesPaths.GENETIC_POPULATION_FILE_NAME
        )

        self.learning_rate_bound = (1, 100)
        self.batch_size_bound = (32, 64)
        self.epochs_bound = (10, 20)
        self.patience_bound = (5, 10)
        self.fine_tune_layers_bound = (1, 3)

        # DEAP stubs for interface compatibility
        self.stats = None
        self.logbook = tools.Logbook()
        self.DEFAULT_LOGBOOK_HEADER = [
            "gen", "nevals", "avg", "std", "min", "max",
            "genotypic_diversity", "phenotypic_diversity"
        ]
        self.logbook.header = self.DEFAULT_LOGBOOK_HEADER
        self.hall_of_fame = tools.HallOfFame(1)

        self.cloud_node = None
        self.evolution_system_metrics = []
        self.evolution_system_metrics_buffer = []

    # ------------------------------------------------------------------
    # Configuration
    # ------------------------------------------------------------------

    def configure_training_parameters_bounds(
        self, lr_min, lr_max, bs_min, bs_max,
        ep_min, ep_max, pa_min, pa_max, ftl_min, ftl_max
    ):
        self.learning_rate_bound = (lr_min, lr_max)
        self.batch_size_bound = (bs_min, bs_max)
        self.epochs_bound = (ep_min, ep_max)
        self.patience_bound = (pa_min, pa_max)
        self.fine_tune_layers_bound = (ftl_min, ftl_max)
        logger.info(
            "Bounds updated: lr=%s bs=%s ep=%s pa=%s ftl=%s",
            self.learning_rate_bound, self.batch_size_bound,
            self.epochs_bound, self.patience_bound, self.fine_tune_layers_bound
        )

    def get_current_training_parameter_bounds(self):
        return {
            "learning_rate_lower_bound": self.learning_rate_bound[0],
            "learning_rate_upper_bound": self.learning_rate_bound[1],
            "batch_size_lower_bound": self.batch_size_bound[0],
            "batch_size_upper_bound": self.batch_size_bound[1],
            "epochs_lower_bound": self.epochs_bound[0],
            "epochs_upper_bound": self.epochs_bound[1],
            "patience_lower_bound": self.patience_bound[0],
            "patience_upper_bound": self.patience_bound[1],
            "fine_tune_layers_lower_bound": self.fine_tune_layers_bound[0],
            "fine_tune_layers_upper_bound": self.fine_tune_layers_bound[1],
        }

    def set_genetic_engine_parameters(self, population_size, number_of_generations, stagnation_limit):
        if population_size is not None:
            self.population_size = population_size
        if number_of_generations is not None:
            self.number_of_generations = number_of_generations
        if stagnation_limit is not None:
            self.stagnation_limit = stagnation_limit
        self.adjust_population_size()

    def get_genetic_engine_parameters(self):
        return {
            "population_size": self.population_size,
            "number_of_generations": self.number_of_generations,
            "stagnation_limit": self.stagnation_limit,
        }

    def set_cloud_node(self, cloud_node):
        self.cloud_node = cloud_node

    # ------------------------------------------------------------------
    # Population management
    # ------------------------------------------------------------------

    def _random_individual(self) -> Individual:
        return Individual.random_individual(
            self.learning_rate_bound,
            self.batch_size_bound,
            self.epochs_bound,
            self.patience_bound,
            self.fine_tune_layers_bound,
        )

    def adjust_population_size(self):
        if self.current_population is None:
            if os.path.exists(self.genetic_population_file_path):
                self.load_population_from_json()
            else:
                return
        current_size = len(self.current_population)
        if self.population_size < current_size:
            self.current_population.sort(
                key=lambda ind: ind.fitness.values[0] if ind.fitness.values else float("inf")
            )
            self.current_population = self.current_population[:self.population_size]
            logger.info("Population truncated to %d.", self.population_size)
        elif self.population_size > current_size:
            num_to_add = self.population_size - current_size
            for _ in range(num_to_add):
                self.current_population.append(self._random_individual())
            logger.info("Added %d individuals; total=%d.", num_to_add, self.population_size)
        else:
            logger.info("Population size unchanged (%d).", self.population_size)

    # ------------------------------------------------------------------
    # Cloud temperature
    # ------------------------------------------------------------------

    def get_cloud_temperature(self):
        if self.cloud_node is None:
            return random.uniform(1, 100)
        url = f"http://{self.cloud_node.ip_address}:{self.cloud_node.port}/cloud/get-cloud-temperature"
        try:
            r = requests.post(url, timeout=None)
            if 200 <= r.status_code < 300:
                temp = r.json().get("cloud_temperature")
                if temp is not None:
                    return temp
        except Exception as e:
            logger.warning("Cloud temperature request failed: %s", e)
        return random.uniform(1, 100)

    # ------------------------------------------------------------------
    # Date helpers
    # ------------------------------------------------------------------

    def set_operating_data_date(self, dates):
        self.operating_data_date = dates

    def clear_operating_data_date(self):
        self.operating_data_date = []

    # ------------------------------------------------------------------
    # Interface stubs
    # ------------------------------------------------------------------

    def should_skip_update(self, fitness, cloud_temperature):
        return False

    def fitness_function(self, individual):
        raise NotImplementedError(
            "fitness_function must not be called on DeepEngine. "
            "Evaluation is handled by the DEEP binary via deep_wrapper_logger.py."
        )

    def setup(self, fog_rabbitmq_host, fog_edge_send_exchange, edge_to_fog_queue):
        self.fog_rabbitmq_host = fog_rabbitmq_host
        self.fog_edge_send_exchange = fog_edge_send_exchange
        self.edge_to_fog_queue = edge_to_fog_queue
        self.logbook.header = self.DEFAULT_LOGBOOK_HEADER
        logger.info("DeepEngine setup complete.")

    # ------------------------------------------------------------------
    # Main evolution entry point
    # ------------------------------------------------------------------

    def evolve(self):
        """
        Run one DEEP optimisation cycle.

        Node selection mirrors GeneticEngine:
          - select_evaluation_node() picks one edge node per round.
          - After the run its timestamp is updated so the next round picks
            a different node, ensuring both edges are used in alternation.

        Environment variables are passed via subprocess.run(env=...) because
        [default_env] is ignored by this DEEP binary build.
        """
        # Reset per-run state
        self.logbook = tools.Logbook()
        self.logbook.header = self.DEFAULT_LOGBOOK_HEADER
        self.stagnation_counter = 0
        self.best_fitness = float("inf")
        self.evolution_system_metrics = []
        self.evolution_system_metrics_buffer = []

        logger.info(
            "Starting DEEP evolution: population_size=%d max_generations=%d",
            self.population_size, self.number_of_generations
        )

        # Clear CSV log from previous run
        self._clear_csv_log()

        # Load previous population (warm-start context)
        if self.current_population is None:
            self.load_population_from_json()
        logger.info(
            "Population before DEEP: %d individuals",
            len(self.current_population) if self.current_population else 0
        )

        # Select evaluation node
        selected_node = select_evaluation_node()
        if selected_node is None:
            logger.error("No evaluation node available; falling back.")
            self._fallback_random()
            return

        validation_url = f"http://{selected_node.ip_address}:{selected_node.port}"
        logger.info("Evaluation node: %s (%s)", selected_node.name, validation_url)

        # Validate prerequisites
        logger.info(
            "Model file: %s exists=%s",
            self.fog_model_file_path,
            os.path.exists(self.fog_model_file_path)
        )
        if not os.path.exists(self.fog_model_file_path):
            logger.error("Model file not found; falling back.")
            self._fallback_random()
            return

        if not os.path.exists(DEEP_BINARY_PATH):
            logger.error("DEEP binary not found at %s; falling back.", DEEP_BINARY_PATH)
            self._fallback_random()
            return

        # Derive date strings
        start_date = (
            self.operating_data_date[0] if len(self.operating_data_date) == 2 else ""
        )
        current_date = (
            self.operating_data_date[1] if len(self.operating_data_date) == 2
            else (self.operating_data_date[0] if len(self.operating_data_date) == 1 else "")
        )

        # Generate DEEP config (bounds + algorithm settings only, no env section)
        config_path = self._generate_deep_config()
        with open(config_path, "r") as f:
            logger.info("DEEP config:\n%s", f.read())

        # Inject wrapper variables via subprocess env.
        # [default_env] in the config is ignored by this DEEP build.
        env = os.environ.copy()
        env["MODEL_FILE_PATH"] = self.fog_model_file_path
        env["VALIDATION_NODE_URL"] = validation_url
        env["CHILD_ID"] = selected_node.id
        env["START_DATE"] = start_date
        env["CURRENT_DATE"] = current_date

        logger.info(
            "Subprocess env: MODEL_FILE_PATH=%s VALIDATION_NODE_URL=%s CHILD_ID=%s "
            "START_DATE=%s CURRENT_DATE=%s",
            env["MODEL_FILE_PATH"], env["VALIDATION_NODE_URL"], env["CHILD_ID"],
            env["START_DATE"], env["CURRENT_DATE"]
        )

        # Launch DEEP
        try:
            logger.info("Launching DEEP: %s --default-name %s", DEEP_BINARY_PATH, config_path)
            result = subprocess.run(
                [DEEP_BINARY_PATH, "--default-name", config_path],
                capture_output=True,
                text=True,
                timeout=3600,
                env=env,
            )
            logger.info("DEEP exit code: %d", result.returncode)
            logger.info("DEEP stdout:\n%s", result.stdout)
            if result.stderr:
                logger.info("DEEP stderr:\n%s", result.stderr)

            if result.returncode != 0:
                logger.error("DEEP non-zero exit; falling back.")
                self._fallback_random()
                return

        except subprocess.TimeoutExpired:
            logger.error("DEEP timed out; falling back.")
            self._fallback_random()
            return
        except Exception as e:
            logger.error("Error running DEEP: %s", e)
            self._fallback_random()
            return
        finally:
            if os.path.exists(config_path):
                os.remove(config_path)

        # Build population from CSV
        best_individuals = self._get_top_n_from_csv(self.population_size)
        if not best_individuals:
            logger.warning("CSV empty after DEEP run; falling back.")
            self._fallback_random()
            return

        self.current_population = best_individuals

        # Pad with random individuals if fewer CSV rows than population_size
        while len(self.current_population) < self.population_size:
            ind = self._random_individual()
            ind.fitness.values = (1e6, 0.0)
            self.current_population.append(ind)

        # Sort best-first, update hall_of_fame (persistent across rounds)
        self.current_population.sort(key=lambda ind: ind.fitness.values[0])
        self.hall_of_fame.update(self.current_population)

        new_best = self.current_population[0].fitness.values[0]
        if new_best < self.best_fitness:
            self.best_fitness = new_best
            self.stagnation_counter = 0
        else:
            self.stagnation_counter += 1

        logger.info(
            "Evolution complete. best_score=%.6f stagnation=%d",
            self.best_fitness, self.stagnation_counter
        )
        logger.info(
            "Hall of fame: %s fitness=%s",
            list(self.hall_of_fame[0]), self.hall_of_fame[0].fitness.values
        )

        self.save_population_to_json()
        self._record_logbook_stub()

        # Update node timestamp — next round selects a different node
        selected_node.last_time_fitness_evaluation_performed_timestamp = time.time_ns()

    # ------------------------------------------------------------------
    # DEEP config generation
    # ------------------------------------------------------------------

    def _generate_deep_config(self) -> str:
        """
        Write a DEEP INI config and return its path.

        Contains only algorithm settings and model bounds.
        Node-specific variables are NOT included — this DEEP build ignores
        [default_env]; they are passed via subprocess env in evolve().

        Verified working format:
          - type=command  (not type=external)
          - Semicolon delimiters on all list fields in [default_model]
          - es_lambda=2  (es_lambda >= 3 caused GLib async-queue crash)
        """
        fd, config_path = tempfile.mkstemp(
            suffix=".conf", dir=SharedResourcesPaths.CACHE_FOLDER_PATH
        )
        os.close(fd)

        lb = ";".join([
            str(self.learning_rate_bound[0]),
            str(self.batch_size_bound[0]),
            str(self.epochs_bound[0]),
            str(self.patience_bound[0]),
            str(self.fine_tune_layers_bound[0]),
        ]) + ";"

        ub = ";".join([
            str(self.learning_rate_bound[1]),
            str(self.batch_size_bound[1]),
            str(self.epochs_bound[1]),
            str(self.patience_bound[1]),
            str(self.fine_tune_layers_bound[1]),
        ]) + ";"

        pv = [
            (self.learning_rate_bound[0] + self.learning_rate_bound[1]) // 2,
            (self.batch_size_bound[0] + self.batch_size_bound[1]) // 2,
            (self.epochs_bound[0] + self.epochs_bound[1]) // 2,
            (self.patience_bound[0] + self.patience_bound[1]) // 2,
            (self.fine_tune_layers_bound[0] + self.fine_tune_layers_bound[1]) // 2,
        ]
        parmvals = ";".join(map(str, pv)) + ";"

        config_content = f"""[default_settings]
max_threads=1
population_size={self.population_size}
max_generations={self.number_of_generations}
es_lambda=3
seed=341640
absolute_iter=10
proportional_stop=1e-04
stop_count=10
logdepth=2
recombination_strategy=de_3_bin
transform=tanh
run_before=initcancel;optpost;optposteval;
run=gcadeep;5;gacdeep;1;gdeep;1;edeep;1;sdeep;24;selde;1;dpupdate;1;checkcancel;1;pdeep;1;substitute;2;optpost;1;optposteval;1;printlog;1;
run_after=optpost;optposteval;printlog;writelog;

[default_target]
debug=0
ignore_cost=0
constrain_aggr=sum
penalty_aggr=sum
mainfunc=target;objfunc;0;1;1;

[default_model]
type=command
partnames=x;y;z;k;l;
partsizes=1;1;1;1;1;
parmvals={parmvals}
parmtypes=1;1;1;1;1;
maskvals=0;0;0;0;0;
tweakvals=1;1;1;1;1;
dparmvals=0;0;0;0;0;
lboundvals={lb}
hboundvals={ub}
limitvals=1;1;1;1;1;
scalevals=1;1;1;1;1;
command=/usr/local/bin/python3 {DEEP_WRAPPER_PATH}
delimiters=\\n=
keys=0;
mapping=0;
num_threads=1
"""
        with open(config_path, "w") as f:
            f.write(config_content)
        logger.info("DEEP config written to %s", config_path)
        return config_path

    # ------------------------------------------------------------------
    # CSV helpers
    # ------------------------------------------------------------------

    def _clear_csv_log(self):
        if os.path.exists(DEEP_CSV_LOG_PATH):
            os.remove(DEEP_CSV_LOG_PATH)
            logger.info("Cleared DEEP CSV log: %s", DEEP_CSV_LOG_PATH)

    def _get_top_n_from_csv(self, n: int) -> List[Individual]:
        """
        Read the CSV written by deep_wrapper_logger and return top-n
        individuals by ascending MSE.  Penalty rows (mse >= 1e5) are skipped.
        """
        if not os.path.exists(DEEP_CSV_LOG_PATH):
            logger.warning("CSV log not found: %s", DEEP_CSV_LOG_PATH)
            return []

        rows = []
        try:
            with open(DEEP_CSV_LOG_PATH, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        mse = float(row["mse"])
                        if not math.isfinite(mse) or mse >= 1e5:
                            continue
                        params = (
                            int(row["lr_raw"]),
                            int(row["batch_size"]),
                            int(row["epochs"]),
                            int(row["patience"]),
                            int(row["fine_tune_layers"]),
                        )
                        rows.append((mse, params))
                    except (KeyError, ValueError) as e:
                        logger.warning("Skipping malformed CSV row %s: %s", row, e)
        except Exception as e:
            logger.error("Error reading CSV log: %s", e)
            return []

        if not rows:
            logger.warning("No valid rows in CSV log.")
            return []

        rows.sort(key=lambda x: x[0])
        individuals = []
        for mse, params in rows[:n]:
            ind = Individual(*params)
            ind.fitness.values = (mse, 0.0)
            individuals.append(ind)

        logger.info(
            "Loaded %d individuals from CSV. best_score=%.6f worst_mse=%.6f",
            len(individuals),
            individuals[0].fitness.values[0],
            individuals[-1].fitness.values[0],
        )
        return individuals

    # ------------------------------------------------------------------
    # Fallback
    # ------------------------------------------------------------------

    def _fallback_random(self):
        logger.warning("Fallback: random population of size %d.", self.population_size)
        self.current_population = [self._random_individual() for _ in range(self.population_size)]
        for ind in self.current_population:
            ind.fitness.values = (1e6, 0.0)
        self.hall_of_fame.update(self.current_population)
        self.save_population_to_json()
        self._record_logbook_stub()

    # ------------------------------------------------------------------
    # Logbook stub
    # ------------------------------------------------------------------

    def _record_logbook_stub(self):
        best_mse = (
            self.current_population[0].fitness.values[0]
            if self.current_population and self.current_population[0].fitness.values
            else 1e6
        )
        self.logbook.record(
            gen=0,
            nevals=len(self.current_population) if self.current_population else 0,
            avg=best_mse,
            std=0.0,
            min=best_mse,
            max=best_mse,
            genotypic_diversity=self.compute_genotypic_diversity(),
            phenotypic_diversity=self.compute_phenotypic_diversity(),
        )

    # ------------------------------------------------------------------
    # Population query
    # ------------------------------------------------------------------

    def get_top_k_individuals(self, k: int) -> List[Individual]:
        if self.current_population is None or len(self.current_population) == 0:
            raise ValueError(
                "Population is empty. Call evolve() before get_top_k_individuals()."
            )
        if k > len(self.current_population):
            logger.warning(
                "Requested %d but population has %d; returning all.",
                k, len(self.current_population)
            )
            k = len(self.current_population)
        return self.current_population[:k]

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save_population_to_json(self):
        if self.current_population is None:
            logger.warning("Nothing to save — population is None.")
            return

        os.makedirs(os.path.dirname(self.genetic_population_file_path), exist_ok=True)

        if os.path.exists(self.genetic_population_file_path):
            os.remove(self.genetic_population_file_path)

        population_data = [
            {
                "chromosome": list(ind),
                "fitness": list(ind.fitness.values) if ind.fitness.values else None,
            }
            for ind in self.current_population
        ]
        header = self.logbook.header or self.DEFAULT_LOGBOOK_HEADER
        logbook_records = [
            [record.get(key) for key in header]
            for record in self.logbook
        ]

        with open(self.genetic_population_file_path, "w") as f:
            json.dump(
                {
                    "population": population_data,
                    "logbook": {"header": header, "records": logbook_records},
                },
                f, indent=2
            )
        logger.info(
            "Saved %d individuals to %s",
            len(self.current_population), self.genetic_population_file_path
        )

    def load_population_from_json(self):
        if not os.path.exists(self.genetic_population_file_path):
            logger.info(
                "No saved population at %s; initialising random population of size %d.",
                self.genetic_population_file_path, self.population_size
            )
            self.current_population = [self._random_individual() for _ in range(self.population_size)]
            return

        logger.info("Loading population from %s", self.genetic_population_file_path)
        with open(self.genetic_population_file_path, "r") as f:
            data = json.load(f)

        population = []
        for ind_data in data.get("population", []):
            chrom = ind_data.get("chromosome", [])
            if len(chrom) == 5:
                ind = Individual(*chrom)
                fitness = ind_data.get("fitness")
                if fitness is not None:
                    ind.fitness.values = tuple(fitness)
                population.append(ind)
            else:
                logger.warning("Skipping invalid chromosome: %s", chrom)

        self.current_population = population if population else [
            self._random_individual() for _ in range(self.population_size)
        ]
        self.hall_of_fame.update(self.current_population)
        logger.info("Loaded %d individuals.", len(self.current_population))

        logbook_data = data.get("logbook")
        if logbook_data and isinstance(logbook_data, dict):
            header = logbook_data.get("header", [])
            records = logbook_data.get("records", [])
            self.logbook = tools.Logbook()
            self.logbook.header = header
            for record in records:
                if isinstance(record, dict):
                    self.logbook.record(**record)
                elif isinstance(record, list):
                    self.logbook.record(**dict(zip(header, record)))
        else:
            self.logbook = tools.Logbook()
            self.logbook.header = self.DEFAULT_LOGBOOK_HEADER

    # ------------------------------------------------------------------
    # Diversity metrics
    # ------------------------------------------------------------------

    def compute_genotypic_diversity(self) -> float:
        if not self.current_population or len(self.current_population) < 2:
            return 0.0
        distances = [
            math.sqrt(sum((a - b) ** 2 for a, b in zip(
                self.current_population[i], self.current_population[j]
            )))
            for i in range(len(self.current_population))
            for j in range(i + 1, len(self.current_population))
        ]
        return sum(distances) / len(distances) if distances else 0.0

    def compute_phenotypic_diversity(self) -> float:
        vals = [
            ind.fitness.values[0]
            for ind in self.current_population
            if ind.fitness.values and math.isfinite(ind.fitness.values[0])
        ]
        if not vals:
            return 0.0
        mean = sum(vals) / len(vals)
        return math.sqrt(sum((x - mean) ** 2 for x in vals) / len(vals))

    def safe_fitness(self, ind) -> float:
        if not ind.fitness.values or not math.isfinite(ind.fitness.values[0]):
            return 1e6
        return sum(ind.fitness.values)