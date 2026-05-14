#!/usr/bin/env python3
"""
deep_wrapper_logger.py
======================
Called by the DEEP binary for each individual evaluation.

DEEP passes 5 hyperparameter values as positional arguments:
    argv[1]  lr_raw            integer (1..100)
    argv[2]  batch_size        integer (32..64)
    argv[3]  epochs            integer (10..20)
    argv[4]  patience          integer (5..10)
    argv[5]  fine_tune_layers  integer (1..3)

Environment variables set by DeepEngine.evolve() via subprocess env:
    MODEL_FILE_PATH      -- absolute path to fog_model.keras
    VALIDATION_NODE_URL  -- base URL of the edge node, e.g. http://edge_app_1:8081
    CHILD_ID             -- id passed as "child_id" in the evaluation payload
    START_DATE           -- optional; start of the operating data window
    CURRENT_DATE         -- current / end date of the operating data window

Output:
    Single float (weighted score) printed to stdout — the only channel DEEP reads.
    One row appended to /app/cache/deep_params_log.csv for post-run ranking.
    Column is named "mse" for backward compatibility with DeepEngine._get_top_n_from_csv().

Fitness logic mirrors GeneticEngine.fitness_function exactly:
    score = min(weighted_score(before_training), weighted_score(after_training))
    weighted_score = sum(weight * metric_value for each metric)

Weights are copied verbatim from shared/utils.py to avoid ImportError
when DEEP spawns this script as an isolated subprocess.
"""

import sys
import os
import base64
import csv
import math
from datetime import datetime

import requests


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PENALTY = 1e6
LOG_FILE = "/app/cache/deep_params_log.csv"
EVAL_ENDPOINT = "/edge/execute-model-evaluation"
REQUEST_TIMEOUT = 300  # seconds

# Bounds — kept in sync with DeepEngine._generate_deep_config() lboundvals/hboundvals
LR_MIN,  LR_MAX  = 1,   100
BS_MIN,  BS_MAX  = 32,   64
EP_MIN,  EP_MAX  = 10,   20
PA_MIN,  PA_MAX  = 5,    10
FTL_MIN, FTL_MAX = 1,     3

# Metric weights — copied verbatim from shared/utils.py.
# Do NOT import shared.utils here: this script runs as an isolated subprocess
# spawned by the DEEP binary and PYTHONPATH may not include /app.
# Keep these values in sync with shared/utils.py manually.
#
#   mse:     0.10
#   logcosh: 0.15  (robust to outliers)
#   huber:   0.10  (less sensitive to outliers than MSE)
#   msle:    0.05  (relative differences)
#   mae:     0.20  (interpretable)
#   r2:     -0.35  (negative: higher R² is better, subtracting lowers the score)
METRIC_WEIGHTS = {
    "mse":     0.10,
    "logcosh": 0.15,
    "huber":   0.10,
    "msle":    0.05,
    "mae":     0.20,
    "r2":     -0.35,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _clamp(value: int, lo: int, hi: int) -> int:
    return max(lo, min(hi, value))


def _compute_weighted_score(metrics: dict) -> float:
    """
    Mirrors GeneticEngine.fitness_function's compute_weighted_score.
    Uses .get(metric, 0.0) to survive missing keys in the edge response.
    """
    score = 0.0
    for metric, weight in METRIC_WEIGHTS.items():
        val = metrics.get(metric)
        if val is None:
            # metric absent — skip rather than crash
            print(f"[wrapper] WARNING: metric '{metric}' missing from response, skipping.", file=sys.stderr)
            continue
        score += weight * float(val)
    return score


def _log(lr_raw, batch_size, epochs, patience, fine_tune_layers, score: float) -> None:
    """
    Append one evaluation row to the CSV log.
    Column header is 'mse' for compatibility with DeepEngine._get_top_n_from_csv().
    The value stored is actually the weighted score — this is intentional.
    """
    file_exists = os.path.isfile(LOG_FILE)
    try:
        os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
        with open(LOG_FILE, "a", newline="") as f:
            writer = csv.writer(f)
            if not file_exists:
                writer.writerow([
                    "timestamp", "lr_raw", "batch_size", "epochs",
                    "patience", "fine_tune_layers", "mse",
                ])
            writer.writerow([
                datetime.now().isoformat(),
                lr_raw, batch_size, epochs, patience, fine_tune_layers,
                score,
            ])
    except Exception as exc:
        print(f"[wrapper] WARNING: could not write CSV log: {exc}", file=sys.stderr)


def _exit_penalty(reason: str = "") -> None:
    """Print penalty to stdout (DEEP reads this) and exit cleanly."""
    if reason:
        print(f"[wrapper] {reason}", file=sys.stderr)
    print(PENALTY, flush=True)
    sys.exit(0)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # 1. Parse command-line arguments from DEEP
    if len(sys.argv) != 6:
        print(f"[wrapper] Expected 5 arguments, got {len(sys.argv) - 1}", file=sys.stderr)
        print(PENALTY, flush=True)
        sys.exit(1)

    try:
        lr_raw           = int(sys.argv[1])
        batch_size       = int(sys.argv[2])
        epochs           = int(sys.argv[3])
        patience         = int(sys.argv[4])
        fine_tune_layers = int(sys.argv[5])
    except ValueError as e:
        _exit_penalty(f"Could not parse arguments: {e}")
        return

    # 2. Clamp to declared bounds (DEEP may slightly overshoot the integer grid)
    lr_raw           = _clamp(lr_raw,           LR_MIN,  LR_MAX)
    batch_size       = _clamp(batch_size,        BS_MIN,  BS_MAX)
    epochs           = _clamp(epochs,            EP_MIN,  EP_MAX)
    patience         = _clamp(patience,          PA_MIN,  PA_MAX)
    fine_tune_layers = _clamp(fine_tune_layers,  FTL_MIN, FTL_MAX)

    learning_rate = lr_raw / 10000.0

    print(
        f"[wrapper] lr={learning_rate} bs={batch_size} ep={epochs} "
        f"pa={patience} ftl={fine_tune_layers}",
        file=sys.stderr,
    )

    # 3. Read environment variables (set by DeepEngine.evolve via subprocess env)
    model_path     = os.environ.get("MODEL_FILE_PATH")
    validation_url = os.environ.get("VALIDATION_NODE_URL")
    child_id       = os.environ.get("CHILD_ID", "deep_wrapper")
    start_date     = os.environ.get("START_DATE")
    current_date   = os.environ.get("CURRENT_DATE")

    if not model_path:
        _exit_penalty("MODEL_FILE_PATH not set.")
        return
    if not validation_url:
        _exit_penalty("VALIDATION_NODE_URL not set.")
        return

    # 4. Read and base64-encode the fog model file
    try:
        with open(model_path, "rb") as f:
            model_b64 = base64.b64encode(f.read()).decode("utf-8")
    except FileNotFoundError:
        _exit_penalty(f"Model file not found: {model_path}")
        return
    except Exception as e:
        _exit_penalty(f"Failed to read model file: {e}")
        return

    # 5. Build evaluation payload — matches GeneticEngine.fitness_function exactly
    payload = {
        "genetic_evaluation": True,
        "start_date":         start_date,
        "current_date":       current_date,
        "is_cache_active":    False,
        "model_type":         "lstm",
        "model_file":         model_b64,
        "learning_rate":      learning_rate,
        "batch_size":         batch_size,
        "epochs":             epochs,
        "patience":           patience,
        "fine_tune_layers":   fine_tune_layers,
        "child_id":           child_id,
        "scope":              1,  # MessageScope.TRAINING
    }

    url = f"{validation_url.rstrip('/')}{EVAL_ENDPOINT}"
    print(f"[wrapper] POST {url}", file=sys.stderr)

    # 6. Send request to edge node
    try:
        resp = requests.post(url, json=payload, timeout=REQUEST_TIMEOUT)
    except requests.exceptions.Timeout:
        _exit_penalty(f"Request to {url} timed out after {REQUEST_TIMEOUT}s.")
        return
    except Exception as e:
        _exit_penalty(f"Request failed: {e}")
        return

    if resp.status_code != 200:
        _exit_penalty(f"Edge returned HTTP {resp.status_code}: {resp.text[:200]}")
        return

    # 7. Parse response
    try:
        data = resp.json()
    except Exception as e:
        _exit_penalty(f"Failed to parse JSON response: {e}")
        return

    # 8. Compute weighted score — mirrors GeneticEngine logic exactly:
    #    score = min(weighted_score(before), weighted_score(after))
    try:
        metrics = data.get("metrics", {})
        before  = metrics.get("before_training", {})
        after   = metrics.get("after_training", {})

        if not before and not after:
            _exit_penalty(f"No metrics in response. Keys: {list(data.keys())}")
            return

        score_before = _compute_weighted_score(before)
        score_after  = _compute_weighted_score(after)
        weighted_score = min(score_before, score_after)

    except Exception as e:
        _exit_penalty(f"Failed to compute weighted score: {e}")
        return

    # 9. Guard against NaN / Inf — DEEP C core may segfault on non-finite floats
    if not math.isfinite(weighted_score):
        _exit_penalty(f"Weighted score is not finite: {weighted_score}")
        return

    print(
        f"[wrapper] score_before={score_before:.6f} score_after={score_after:.6f} "
        f"weighted_score={weighted_score:.6f}",
        file=sys.stderr,
    )

    # 10. Log to CSV for DeepEngine._get_top_n_from_csv()
    _log(lr_raw, batch_size, epochs, patience, fine_tune_layers, weighted_score)

    # 11. Return score to DEEP (minimised)
    print(weighted_score, flush=True)
    sys.exit(0)


if __name__ == "__main__":
    main()