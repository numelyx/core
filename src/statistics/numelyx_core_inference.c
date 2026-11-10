/**
 * @file numelyx_core_inference.c
 * @brief Delta-method and Wald tests using stable covariance factorizations.
 */
#include "numelyx_core_inference.h"
#include "numelyx_core_size.h"
#include "numelyx_core_distributions.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

static void inference_kahan_add(double *sum, double *correction, double value)
{
    double adjusted = value - *correction;
    double next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

numelyx_core_inference_status numelyx_core_delta_method_variance(
    const double *gradient, const double *covariance,
    size_t parameter_count, double *variance
)
{
    size_t matrix_elements;
    size_t i;
    size_t j;
    double sum = 0.0;
    double correction = 0.0;
    double absolute_sum = 0.0;

    if (gradient == NULL || covariance == NULL || variance == NULL ||
        parameter_count == 0)
        return NUMELYX_INFERENCE_INVALID_ARGUMENT;
    if (!numelyx_core_size_multiply(
            parameter_count, parameter_count, &matrix_elements))
        return NUMELYX_INFERENCE_SIZE_OVERFLOW;

    for (i = 0; i < parameter_count; i++) {
        double row = 0.0;
        if (!isfinite(gradient[i]))
            return NUMELYX_INFERENCE_NON_FINITE_INPUT;
        for (j = 0; j < parameter_count; j++) {
            double value = covariance[i * parameter_count + j];
            if (!isfinite(value))
                return NUMELYX_INFERENCE_NON_FINITE_INPUT;
            row += value * gradient[j];
            if (!isfinite(row))
                return NUMELYX_INFERENCE_NUMERICAL_FAILURE;
        }
        {
            double term = gradient[i] * row;
            inference_kahan_add(&sum, &correction, term);
            absolute_sum += fabs(term);
            if (!isfinite(sum) || !isfinite(absolute_sum))
                return NUMELYX_INFERENCE_NUMERICAL_FAILURE;
        }
    }
    if (sum < 0.0) {
        double tolerance = 64.0 * DBL_EPSILON * fmax(1.0, absolute_sum);
        if (sum < -tolerance)
            return NUMELYX_INFERENCE_NUMERICAL_FAILURE;
        sum = 0.0;
    }
    *variance = sum;
    return NUMELYX_INFERENCE_SUCCESS;
}

numelyx_core_inference_status numelyx_core_delta_method_standard_error(
    const double *gradient, const double *covariance,
    size_t parameter_count, double *standard_error
)
{
    double variance;
    numelyx_core_inference_status status;
    if (standard_error == NULL) return NUMELYX_INFERENCE_INVALID_ARGUMENT;
    status = numelyx_core_delta_method_variance(
        gradient, covariance, parameter_count, &variance);
    if (status != NUMELYX_INFERENCE_SUCCESS) return status;
    *standard_error = sqrt(variance);
    return NUMELYX_INFERENCE_SUCCESS;
}

static bool cholesky_factor(
    const double *matrix, double *factor, size_t count
)
{
    size_t i;
    size_t j;
    size_t k;
    double scale = 0.0;
    double tolerance;

    for (i = 0; i < count; i++) {
        if (matrix[i * count + i] > scale)
            scale = matrix[i * count + i];
    }
    if (!isfinite(scale) || scale <= 0.0) return false;
    tolerance = 64.0 * DBL_EPSILON * scale;
    memset(factor, 0, count * count * sizeof(*factor));

    for (i = 0; i < count; i++) {
        for (j = 0; j <= i; j++) {
            double value = matrix[i * count + j];
            for (k = 0; k < j; k++)
                value -= factor[i * count + k] * factor[j * count + k];
            if (!isfinite(value)) return false;
            if (i == j) {
                if (value <= tolerance) return false;
                factor[i * count + j] = sqrt(value);
            } else {
                factor[i * count + j] =
                    value / factor[j * count + j];
                if (!isfinite(factor[i * count + j])) return false;
            }
        }
    }
    return true;
}

numelyx_core_inference_status numelyx_core_wald_test(
    const double *estimate, const double *null_value,
    const double *restriction_matrix, const double *covariance,
    size_t parameter_count, size_t restriction_count,
    double *wald_statistic, double *log_p_value
)
{
    size_t covariance_elements;
    size_t restriction_elements;
    size_t restriction_square;
    size_t covariance_bytes;
    size_t restriction_covariance_bytes;
    size_t restriction_bytes;
    size_t vector_bytes;
    size_t i;
    size_t j;
    size_t k;
    double *restricted_covariance = NULL;
    double *factor = NULL;
    double *rv = NULL;
    double *difference = NULL;
    double *forward = NULL;
    double *solution = NULL;
    double statistic = 0.0;
    double statistic_correction = 0.0;
    double covariance_scale = 0.0;
    numelyx_core_inference_status status =
        NUMELYX_INFERENCE_NUMERICAL_FAILURE;

    if (estimate == NULL || null_value == NULL ||
        restriction_matrix == NULL || covariance == NULL ||
        wald_statistic == NULL || log_p_value == NULL ||
        parameter_count == 0 || restriction_count == 0)
        return NUMELYX_INFERENCE_INVALID_ARGUMENT;
    if (!numelyx_core_size_multiply(
            parameter_count, parameter_count, &covariance_elements) ||
        !numelyx_core_size_multiply(
            restriction_count, parameter_count, &restriction_elements) ||
        !numelyx_core_size_multiply(
            restriction_count, restriction_count, &restriction_square) ||
        !numelyx_core_size_multiply(
            covariance_elements, sizeof(double), &covariance_bytes) ||
        !numelyx_core_size_multiply(
            restriction_square, sizeof(double), &restriction_covariance_bytes) ||
        !numelyx_core_size_multiply(
            restriction_elements, sizeof(double), &restriction_bytes) ||
        !numelyx_core_size_multiply(
            restriction_count, sizeof(double), &vector_bytes))
        return NUMELYX_INFERENCE_SIZE_OVERFLOW;

    restricted_covariance = (double *) calloc(1, restriction_covariance_bytes);
    factor = (double *) malloc(restriction_covariance_bytes);
    rv = (double *) malloc(restriction_bytes);
    difference = (double *) malloc(vector_bytes);
    forward = (double *) malloc(vector_bytes);
    solution = (double *) malloc(vector_bytes);
    if (restricted_covariance == NULL || factor == NULL || rv == NULL ||
        difference == NULL || forward == NULL || solution == NULL) {
        status = NUMELYX_INFERENCE_OUT_OF_MEMORY;
        goto cleanup;
    }

    for (i = 0; i < parameter_count; i++) {
        if (!isfinite(estimate[i]))
            { status = NUMELYX_INFERENCE_NON_FINITE_INPUT; goto cleanup; }
        for (j = 0; j < parameter_count; j++) {
            double value = covariance[i * parameter_count + j];
            if (!isfinite(value))
                { status = NUMELYX_INFERENCE_NON_FINITE_INPUT; goto cleanup; }
            if (fabs(value) > covariance_scale)
                covariance_scale = fabs(value);
            if (fabs(value - covariance[j * parameter_count + i]) >
                1e-10 * fmax(1.0, covariance_scale)) {
                status = NUMELYX_INFERENCE_INVALID_ARGUMENT;
                goto cleanup;
            }
        }
    }
    for (i = 0; i < restriction_count; i++) {
        double value = -null_value[i];
        double correction = 0.0;
        if (!isfinite(null_value[i])) {
            status = NUMELYX_INFERENCE_NON_FINITE_INPUT;
            goto cleanup;
        }
        for (j = 0; j < parameter_count; j++) {
            double r = restriction_matrix[i * parameter_count + j];
            if (!isfinite(r)) {
                status = NUMELYX_INFERENCE_NON_FINITE_INPUT;
                goto cleanup;
            }
            inference_kahan_add(&value, &correction, r * estimate[j]);
        }
        difference[i] = value;
    }

    for (i = 0; i < restriction_count; i++) {
        for (j = 0; j < parameter_count; j++) {
            double value = 0.0;
            double correction = 0.0;
            for (k = 0; k < parameter_count; k++) {
                inference_kahan_add(&value, &correction,
                    restriction_matrix[i * parameter_count + k] *
                    covariance[k * parameter_count + j]);
            }
            if (!isfinite(value)) {
                status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
                goto cleanup;
            }
            rv[i * parameter_count + j] = value;
        }
    }
    for (i = 0; i < restriction_count; i++) {
        for (j = 0; j < restriction_count; j++) {
            double value = 0.0;
            double correction = 0.0;
            for (k = 0; k < parameter_count; k++) {
                inference_kahan_add(&value, &correction,
                    rv[i * parameter_count + k] *
                    restriction_matrix[j * parameter_count + k]);
            }
            restricted_covariance[i * restriction_count + j] = value;
        }
    }
    for (i = 0; i < restriction_square; i++) {
        if (!isfinite(restricted_covariance[i])) {
            status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
            goto cleanup;
        }
    }
    for (i = 0; i < restriction_count; i++) {
        for (j = i + 1; j < restriction_count; j++) {
            double symmetric =
                0.5 * restricted_covariance[i * restriction_count + j] +
                0.5 * restricted_covariance[j * restriction_count + i];
            restricted_covariance[i * restriction_count + j] = symmetric;
            restricted_covariance[j * restriction_count + i] = symmetric;
        }
    }
    if (!cholesky_factor(
            restricted_covariance, factor, restriction_count)) {
        status = NUMELYX_INFERENCE_SINGULAR_COVARIANCE;
        goto cleanup;
    }

    for (i = 0; i < restriction_count; i++) {
        double value = difference[i];
        for (j = 0; j < i; j++)
            value -= factor[i * restriction_count + j] * forward[j];
        forward[i] = value / factor[i * restriction_count + i];
        if (!isfinite(forward[i])) {
            status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
            goto cleanup;
        }
    }
    for (i = restriction_count; i-- > 0;) {
        double value = forward[i];
        for (j = i + 1; j < restriction_count; j++)
            value -= factor[j * restriction_count + i] * solution[j];
        solution[i] = value / factor[i * restriction_count + i];
        if (!isfinite(solution[i])) {
            status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
            goto cleanup;
        }
    }
    for (i = 0; i < restriction_count; i++) {
        inference_kahan_add(&statistic, &statistic_correction,
            difference[i] * solution[i]);
    }
    if (!isfinite(statistic) || statistic < -1e-10) {
        status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
        goto cleanup;
    }
    if (statistic < 0.0) statistic = 0.0;
    {
        double log_probability;
        numelyx_core_numeric_status numeric_status =
            numelyx_core_chi_square_log_survival(
                statistic, (double) restriction_count, &log_probability);
        if (numeric_status != NUMERUS_NUMERIC_SUCCESS) {
            status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
            goto cleanup;
        }
        *wald_statistic = statistic;
        *log_p_value = log_probability;
    }
    status = NUMELYX_INFERENCE_SUCCESS;

cleanup:
    free(restricted_covariance);
    free(factor);
    free(rv);
    free(difference);
    free(forward);
    free(solution);
    return status;
}
