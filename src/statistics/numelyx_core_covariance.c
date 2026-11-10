/**
 * @file numelyx_core_covariance.c
 * @brief HC, cluster-robust, and Newey-West sandwich covariance estimators.
 */
#include "numelyx_core_covariance.h"
#include "numelyx_core_size.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

static void kahan_add(double *sum, double *correction, double value)
{
    double adjusted = value - *correction;
    double next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

static bool validate_covariance_inputs(
    const double *bread, const double *scores,
    size_t observations, size_t parameters, double *output
)
{
    size_t i;
    size_t j;
    if (bread == NULL || scores == NULL || output == NULL ||
        observations == 0 || parameters == 0) return false;
    for (i = 0; i < parameters * parameters; i++) {
        if (!isfinite(bread[i])) return false;
    }
    for (i = 0; i < observations; i++) {
        for (j = 0; j < parameters; j++) {
            if (!isfinite(scores[i * parameters + j])) return false;
        }
    }
    return true;
}

static numelyx_core_covariance_status finish_sandwich(
    const double *bread, const double *meat, size_t parameters,
    double correction_factor, double *covariance,
    double *temporary, double *candidate
)
{
    size_t i;
    size_t j;
    size_t k;
    for (i = 0; i < parameters; i++) {
        for (j = 0; j < parameters; j++) {
            double value = 0.0;
            for (k = 0; k < parameters; k++) {
                value += bread[i * parameters + k] *
                    meat[k * parameters + j];
                if (!isfinite(value))
                    return NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
            }
            temporary[i * parameters + j] = value * correction_factor;
            if (!isfinite(temporary[i * parameters + j]))
                return NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
        }
    }
    for (i = 0; i < parameters; i++) {
        for (j = 0; j < parameters; j++) {
            double value = 0.0;
            for (k = 0; k < parameters; k++) {
                value += temporary[i * parameters + k] *
                    bread[j * parameters + k];
                if (!isfinite(value))
                    return NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
            }
            candidate[i * parameters + j] = value;
        }
    }
    for (i = 0; i < parameters; i++) {
        for (j = i; j < parameters; j++) {
            double symmetric =
                0.5 * candidate[i * parameters + j] +
                0.5 * candidate[j * parameters + i];
            if (!isfinite(symmetric))
                return NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
            candidate[i * parameters + j] = symmetric;
            candidate[j * parameters + i] = symmetric;
        }
    }
    memcpy(covariance, candidate,
        parameters * parameters * sizeof(*covariance));
    return NUMELYX_COVARIANCE_SUCCESS;
}

static numelyx_core_covariance_status allocate_workspaces(
    size_t parameters, double **meat, double **compensation,
    double **temporary, double **candidate, size_t *matrix_bytes
)
{
    size_t elements;
    if (!numelyx_core_size_multiply(parameters, parameters, &elements) ||
        !numelyx_core_size_multiply(elements, sizeof(double), matrix_bytes))
        return NUMELYX_COVARIANCE_SIZE_OVERFLOW;
    *meat = (double *) calloc(1, *matrix_bytes);
    *compensation = (double *) calloc(1, *matrix_bytes);
    *temporary = (double *) malloc(*matrix_bytes);
    *candidate = (double *) malloc(*matrix_bytes);
    if (*meat == NULL || *compensation == NULL ||
        *temporary == NULL || *candidate == NULL)
        return NUMELYX_COVARIANCE_OUT_OF_MEMORY;
    return NUMELYX_COVARIANCE_SUCCESS;
}

static void free_workspaces(
    double *meat, double *compensation, double *temporary, double *candidate
)
{
    free(meat);
    free(compensation);
    free(temporary);
    free(candidate);
}

numelyx_core_covariance_status numelyx_core_sandwich_covariance(
    const double *bread, const double *scores,
    size_t observations, size_t parameters,
    numelyx_core_covariance_correction correction,
    double *covariance
)
{
    double *meat = NULL;
    double *compensation = NULL;
    double *temporary = NULL;
    double *candidate = NULL;
    size_t matrix_bytes = 0;
    size_t score_elements;
    size_t observation;
    size_t i;
    size_t j;
    double correction_factor = 1.0;
    numelyx_core_covariance_status status;

    if (correction != NUMELYX_COVARIANCE_HC0 &&
        correction != NUMELYX_COVARIANCE_HC1)
        return NUMELYX_COVARIANCE_INVALID_ARGUMENT;
    if (bread == NULL || scores == NULL || covariance == NULL ||
        observations == 0 || parameters == 0)
        return NUMELYX_COVARIANCE_INVALID_ARGUMENT;
    if (correction == NUMELYX_COVARIANCE_HC1) {
        if (observations <= parameters)
            return NUMELYX_COVARIANCE_DOMAIN_ERROR;
        correction_factor = (double) observations /
            (double) (observations - parameters);
    }
    if (!numelyx_core_size_multiply(
            observations, parameters, &score_elements))
        return NUMELYX_COVARIANCE_SIZE_OVERFLOW;
    status = allocate_workspaces(parameters, &meat, &compensation,
        &temporary, &candidate, &matrix_bytes);
    if (status != NUMELYX_COVARIANCE_SUCCESS) goto cleanup;
    if (!validate_covariance_inputs(
            bread, scores, observations, parameters, covariance)) {
        status = NUMELYX_COVARIANCE_NON_FINITE_INPUT;
        goto cleanup;
    }

    for (observation = 0; observation < observations; observation++) {
        const double *score = scores + observation * parameters;
        for (i = 0; i < parameters; i++) {
            for (j = 0; j < parameters; j++) {
                size_t index = i * parameters + j;
                kahan_add(&meat[index], &compensation[index],
                    score[i] * score[j]);
                if (!isfinite(meat[index])) {
                    status = NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
                    goto cleanup;
                }
            }
        }
    }
    status = finish_sandwich(bread, meat, parameters, correction_factor,
        covariance, temporary, candidate);

cleanup:
    free_workspaces(meat, compensation, temporary, candidate);
    return status;
}

numelyx_core_covariance_status numelyx_core_cluster_sandwich_covariance(
    const double *bread, const double *scores, const size_t *cluster_ids,
    size_t observations, size_t parameters, size_t cluster_count,
    numelyx_core_cluster_correction correction, double *covariance
)
{
    double *cluster_scores = NULL;
    double *cluster_compensation = NULL;
    unsigned char *cluster_seen = NULL;
    double *meat = NULL;
    double *compensation = NULL;
    double *temporary = NULL;
    double *candidate = NULL;
    size_t matrix_bytes = 0;
    size_t score_elements;
    size_t cluster_elements;
    size_t cluster_bytes;
    size_t cluster_seen_bytes;
    size_t observation;
    size_t i;
    size_t j;
    double correction_factor = 1.0;
    numelyx_core_covariance_status status;

    if (cluster_ids == NULL || cluster_count == 0 ||
        (correction != NUMELYX_COVARIANCE_CLUSTER_CR0 &&
         correction != NUMELYX_COVARIANCE_CLUSTER_CR1))
        return NUMELYX_COVARIANCE_INVALID_ARGUMENT;
    if (correction == NUMELYX_COVARIANCE_CLUSTER_CR1) {
        if (cluster_count <= 1) return NUMELYX_COVARIANCE_DOMAIN_ERROR;
        correction_factor = (double) cluster_count /
            (double) (cluster_count - 1);
    }
    if (bread == NULL || scores == NULL || covariance == NULL ||
        observations == 0 || parameters == 0 || cluster_count > observations)
        return NUMELYX_COVARIANCE_INVALID_ARGUMENT;
    if (!numelyx_core_size_multiply(
            observations, parameters, &score_elements) ||
        !numelyx_core_size_multiply(
            cluster_count, parameters, &cluster_elements) ||
        !numelyx_core_size_multiply(
            cluster_elements, sizeof(double), &cluster_bytes) ||
        !numelyx_core_size_multiply(
            cluster_count, sizeof(*cluster_seen), &cluster_seen_bytes)) {
        return NUMELYX_COVARIANCE_SIZE_OVERFLOW;
    }
    status = allocate_workspaces(parameters, &meat, &compensation,
        &temporary, &candidate, &matrix_bytes);
    if (status != NUMELYX_COVARIANCE_SUCCESS) goto cleanup;
    cluster_scores = (double *) calloc(1, cluster_bytes);
    cluster_compensation = (double *) calloc(1, cluster_bytes);
    cluster_seen = (unsigned char *) calloc(1, cluster_seen_bytes);
    if (cluster_scores == NULL || cluster_compensation == NULL ||
        cluster_seen == NULL) {
        status = NUMELYX_COVARIANCE_OUT_OF_MEMORY;
        goto cleanup;
    }
    if (!validate_covariance_inputs(
            bread, scores, observations, parameters, covariance)) {
        status = NUMELYX_COVARIANCE_NON_FINITE_INPUT;
        goto cleanup;
    }

    for (observation = 0; observation < observations; observation++) {
        size_t cluster = cluster_ids[observation];
        if (cluster >= cluster_count) {
            status = NUMELYX_COVARIANCE_DOMAIN_ERROR;
            goto cleanup;
        }
        cluster_seen[cluster] = 1;
        for (i = 0; i < parameters; i++) {
            size_t index = cluster * parameters + i;
            kahan_add(&cluster_scores[index], &cluster_compensation[index],
                scores[observation * parameters + i]);
            if (!isfinite(cluster_scores[index])) {
                status = NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
                goto cleanup;
            }
        }
    }
    for (i = 0; i < cluster_count; i++) {
        const double *score;
        if (cluster_seen[i] == 0) {
            status = NUMELYX_COVARIANCE_DOMAIN_ERROR;
            goto cleanup;
        }
        score = cluster_scores + i * parameters;
        for (j = 0; j < parameters; j++) {
            size_t k;
            for (k = 0; k < parameters; k++) {
                size_t index = j * parameters + k;
                kahan_add(&meat[index], &compensation[index],
                    score[j] * score[k]);
                if (!isfinite(meat[index])) {
                    status = NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
                    goto cleanup;
                }
            }
        }
    }
    status = finish_sandwich(bread, meat, parameters, correction_factor,
        covariance, temporary, candidate);

cleanup:
    free(cluster_scores);
    free(cluster_compensation);
    free(cluster_seen);
    free_workspaces(meat, compensation, temporary, candidate);
    return status;
}

numelyx_core_covariance_status numelyx_core_hac_sandwich_covariance(
    const double *bread, const double *scores, size_t observations,
    size_t parameters, size_t bandwidth, double *covariance
)
{
    double *meat = NULL;
    double *compensation = NULL;
    double *temporary = NULL;
    double *candidate = NULL;
    size_t matrix_bytes = 0;
    size_t score_elements;
    size_t observation;
    size_t lag;
    size_t i;
    size_t j;
    numelyx_core_covariance_status status;

    if (bread == NULL || scores == NULL || covariance == NULL ||
        observations == 0 || parameters == 0)
        return NUMELYX_COVARIANCE_INVALID_ARGUMENT;
    if (!numelyx_core_size_multiply(
            observations, parameters, &score_elements))
        return NUMELYX_COVARIANCE_SIZE_OVERFLOW;
    status = allocate_workspaces(parameters, &meat, &compensation,
        &temporary, &candidate, &matrix_bytes);
    if (status != NUMELYX_COVARIANCE_SUCCESS) goto cleanup;
    if (!validate_covariance_inputs(
            bread, scores, observations, parameters, covariance)) {
        status = NUMELYX_COVARIANCE_NON_FINITE_INPUT;
        goto cleanup;
    }
    if (bandwidth >= observations) bandwidth = observations - 1;

    for (observation = 0; observation < observations; observation++) {
        const double *score = scores + observation * parameters;
        for (i = 0; i < parameters; i++) {
            for (j = 0; j < parameters; j++) {
                size_t index = i * parameters + j;
                kahan_add(&meat[index], &compensation[index],
                    score[i] * score[j]);
            }
        }
    }

    for (lag = 1; lag <= bandwidth; lag++) {
        double weight = 1.0 - (double) lag / (double) (bandwidth + 1);
        for (observation = lag; observation < observations; observation++) {
            const double *current = scores + observation * parameters;
            const double *lagged = scores + (observation - lag) * parameters;
            for (i = 0; i < parameters; i++) {
                for (j = 0; j < parameters; j++) {
                    size_t index = i * parameters + j;
                    kahan_add(&meat[index], &compensation[index],
                        weight * current[i] * lagged[j]);
                    kahan_add(&meat[j * parameters + i],
                        &compensation[j * parameters + i],
                        weight * lagged[j] * current[i]);
                }
            }
        }
    }
    for (i = 0; i < parameters * parameters; i++) {
        if (!isfinite(meat[i])) {
            status = NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
            goto cleanup;
        }
    }
    status = finish_sandwich(bread, meat, parameters, 1.0,
        covariance, temporary, candidate);

cleanup:
    free_workspaces(meat, compensation, temporary, candidate);
    return status;
}


numelyx_core_covariance_status numelyx_core_hc_sandwich_covariance(
    const double *bread,
    const double *design,
    const double *residuals,
    const double *leverage,
    size_t observations,
    size_t parameters,
    numelyx_core_hc_type type,
    double *covariance
)
{
    size_t score_elements;
    size_t score_bytes;
    size_t covariance_elements;
    size_t covariance_bytes;
    size_t observation;
    size_t parameter;
    double *scores = NULL;
    double *candidate = NULL;
    numelyx_core_covariance_correction correction =
        NUMELYX_COVARIANCE_HC0;
    numelyx_core_covariance_status status;

    if (bread == NULL || design == NULL || residuals == NULL ||
        covariance == NULL || observations == 0 || parameters == 0)
        return NUMELYX_COVARIANCE_INVALID_ARGUMENT;
    if (type < NUMELYX_COVARIANCE_HC_TYPE_0 ||
        type > NUMELYX_COVARIANCE_HC_TYPE_3)
        return NUMELYX_COVARIANCE_INVALID_ARGUMENT;
    if ((type == NUMELYX_COVARIANCE_HC_TYPE_2 ||
         type == NUMELYX_COVARIANCE_HC_TYPE_3) && leverage == NULL)
        return NUMELYX_COVARIANCE_INVALID_ARGUMENT;
    if (type == NUMELYX_COVARIANCE_HC_TYPE_1)
        correction = NUMELYX_COVARIANCE_HC1;

    if (!numelyx_core_size_multiply(
            observations, parameters, &score_elements) ||
        !numelyx_core_size_multiply(
            score_elements, sizeof(double), &score_bytes) ||
        !numelyx_core_size_multiply(
            parameters, parameters, &covariance_elements) ||
        !numelyx_core_size_multiply(
            covariance_elements, sizeof(double), &covariance_bytes))
        return NUMELYX_COVARIANCE_SIZE_OVERFLOW;

    scores = (double *) malloc(score_bytes);
    candidate = (double *) malloc(covariance_bytes);
    if (scores == NULL || candidate == NULL) {
        status = NUMELYX_COVARIANCE_OUT_OF_MEMORY;
        goto cleanup;
    }

    for (observation = 0; observation < observations; observation++) {
        double adjustment = 1.0;
        if (!isfinite(residuals[observation])) {
            status = NUMELYX_COVARIANCE_NON_FINITE_INPUT;
            goto cleanup;
        }
        if (type == NUMELYX_COVARIANCE_HC_TYPE_2 ||
            type == NUMELYX_COVARIANCE_HC_TYPE_3) {
            double remaining;
            if (!isfinite(leverage[observation])) {
                status = NUMELYX_COVARIANCE_NON_FINITE_INPUT;
                goto cleanup;
            }
            if (leverage[observation] < 0.0 ||
                leverage[observation] >= 1.0) {
                status = NUMELYX_COVARIANCE_DOMAIN_ERROR;
                goto cleanup;
            }
            remaining = 1.0 - leverage[observation];
            adjustment = type == NUMELYX_COVARIANCE_HC_TYPE_2
                ? sqrt(remaining)
                : remaining;
        }
        for (parameter = 0; parameter < parameters; parameter++) {
            double x = design[observation * parameters + parameter];
            double score;
            if (!isfinite(x)) {
                status = NUMELYX_COVARIANCE_NON_FINITE_INPUT;
                goto cleanup;
            }
            score = x * (residuals[observation] / adjustment);
            if (!isfinite(score)) {
                status = NUMELYX_COVARIANCE_NUMERICAL_FAILURE;
                goto cleanup;
            }
            scores[observation * parameters + parameter] = score;
        }
    }

    status = numelyx_core_sandwich_covariance(
        bread, scores, observations, parameters, correction, candidate);
    if (status == NUMELYX_COVARIANCE_SUCCESS) {
        memcpy(covariance, candidate, covariance_bytes);
    }

cleanup:
    free(scores);
    free(candidate);
    return status;
}
