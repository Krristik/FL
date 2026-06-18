import os
import numpy as np
import tensorflow as tf

from cloud_node.cloud_resources_paths import CloudResourcesPaths
from shared.shared_resources_paths import SharedResourcesPaths
from shared.logging_config import logger
from model_architectures import (create_initial_lstm_model, create_attention_lstm_model,
                                 create_enhanced_attention_lstm_model)

# Suppress TensorFlow low-level logs
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Drift detection threshold factor; fog models with drift > (threshold_factor * median) are considered outliers
DRIFT_THRESHOLD_FACTOR = 3


def create_model(model_type):
    if model_type == "1":
        return create_initial_lstm_model()
    elif model_type == "2":
        return create_attention_lstm_model()
    elif model_type == "3":
        return create_enhanced_attention_lstm_model()
    else:
        logger.warning(f"Unknown model_type '{model_type}', defaulting to initial LSTM model.")
        return create_initial_lstm_model()


def aggregate_fog_models(received_fog_messages: dict):
    logger.info("Running model aggregation in cloud node with drift detection and robust correction.")

    custom_objects = {
        "LogCosh": tf.keras.losses.LogCosh(),
        "mse": tf.keras.losses.MeanSquaredError(),
        "Huber": tf.keras.losses.Huber()
    }
    cloud_model_file_path = os.path.join(
        CloudResourcesPaths.MODELS_FOLDER_PATH,
        CloudResourcesPaths.CLOUD_MODEL_FILE_NAME
    )

    cache_cloud_model_file_path = os.path.join(SharedResourcesPaths.CACHE_FOLDER_PATH,
                                               CloudResourcesPaths.CLOUD_MODEL_FILE_NAME)

    if os.path.exists(cloud_model_file_path):
        cloud_model = tf.keras.models.load_model(cloud_model_file_path, custom_objects=custom_objects)
    else:
        cloud_model = tf.keras.models.load_model(cache_cloud_model_file_path, custom_objects=custom_objects)
    cloud_model_weights = cloud_model.get_weights()

    fog_models_data = []
    lambda_values = []    # cumulative performance factors (λ) — для совместимости
    dataset_sizes = []    # N_c = Σ|D_k| для каждого фога

    for fog_id, message in received_fog_messages.items():
        fog_model_path = message.get("fog_model_file_path")
        fog_model = tf.keras.models.load_model(fog_model_path, custom_objects=custom_objects)
        fog_model_weights = fog_model.get_weights()
        lambda_prev_value = float(message.get("lambda_prev"))
        dataset_size = float(message.get("dataset_size", 1))  # N_c, fallback=1
        fog_models_data.append(fog_model_weights)
        lambda_values.append(lambda_prev_value)
        dataset_sizes.append(dataset_size)

    # Normalize lambda values (для совместимости и логирования)
    lambda_sum = sum(lambda_values)
    normalized_lambdas = [lmbd / lambda_sum for lmbd in lambda_values] if lambda_sum != 0 else [0 for _ in lambda_values]

    # μ_c = N_c / Σ N_j  (формула 2.12)
    total_dataset_size = sum(dataset_sizes)
    if total_dataset_size > 0:
        mu_values = [n / total_dataset_size for n in dataset_sizes]
    else:
        logger.warning("Total dataset size is 0, falling back to equal weights.")
        mu_values = [1.0 / len(dataset_sizes) for _ in dataset_sizes]

    logger.info("Lambda participation counters and dataset-based weights received from Fog nodes:")
    for i, fog_id in enumerate(received_fog_messages.keys()):
        logger.info(f"  {fog_id}: lambda_prev={lambda_values[i]}, normalized_beta={normalized_lambdas[i]:.4f}, "
                    f"dataset_size={dataset_sizes[i]}, mu={mu_values[i]:.4f}")

    aggregated_weights = [np.zeros_like(layer) for layer in cloud_model_weights]

    for layer_index in range(len(aggregated_weights)):
        cloud_layer = cloud_model_weights[layer_index]

        # w(t) = β * w(t-1) + (1 - β) * Σ μ_c * w_c(t)  (формула 2.11, β=0.5)
        fog_aggregate = sum(
            mu_values[i] * fog_models_data[i][layer_index]
            for i in range(len(fog_models_data))
        )
        aggregated_layer = cloud_layer * 0.5 + fog_aggregate * 0.5
        aggregated_weights[layer_index] = aggregated_layer

    cloud_model.set_weights(aggregated_weights)
    cloud_model.save(cloud_model_file_path)
    cloud_model.save(cache_cloud_model_file_path)

    logger.info(f"Cloud model aggregation correction completed. New cloud model saved at: "
                f"{cloud_model_file_path} and cached.")
