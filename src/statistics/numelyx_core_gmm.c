/**
 * @file numelyx_core_gmm.c
 * @brief Allocation-checked generic GMM objective and derivative evaluation.
 */
#include "numelyx_core_gmm.h"
#include "numelyx_core_distributions.h"
#include "numelyx_core_size.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static void gmm_kahan_add(double *sum, double *correction, double value)
{
    double adjusted = value - *correction;
    double next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

numelyx_core_inference_status numelyx_core_gmm_evaluate(
    const double *moments,
    const double *weighting_matrix,
    const double *moment_jacobian,
    size_t moment_count,
    size_t parameter_count,
    double *objective,
    double *gradient,
    double *gauss_newton_hessian
)
{
    size_t moment_square;
    size_t jacobian_elements;
    size_t parameter_square;
    size_t moment_bytes;
    size_t jacobian_bytes;
    size_t parameter_bytes;
    size_t parameter_matrix_bytes;
    size_t i;
    size_t j;
    size_t k;
    double *weighted_moments = NULL;
    double *weighted_jacobian = NULL;
    double *gradient_candidate = NULL;
    double *hessian_candidate = NULL;
    double objective_value = 0.0;
    double objective_correction = 0.0;
    double objective_absolute_sum = 0.0;
    double weighting_scale = 0.0;
    numelyx_core_inference_status status =
        NUMELYX_INFERENCE_NUMERICAL_FAILURE;

    if (moments == NULL || weighting_matrix == NULL ||
        moment_jacobian == NULL || objective == NULL ||
        moment_count == 0 || parameter_count == 0)
        return NUMELYX_INFERENCE_INVALID_ARGUMENT;
    if (!numelyx_core_size_multiply(
            moment_count, moment_count, &moment_square) ||
        !numelyx_core_size_multiply(
            moment_count, parameter_count, &jacobian_elements) ||
        !numelyx_core_size_multiply(
            parameter_count, parameter_count, &parameter_square) ||
        !numelyx_core_size_multiply(
            moment_count, sizeof(double), &moment_bytes) ||
        !numelyx_core_size_multiply(
            jacobian_elements, sizeof(double), &jacobian_bytes) ||
        !numelyx_core_size_multiply(
            parameter_count, sizeof(double), &parameter_bytes) ||
        !numelyx_core_size_multiply(
            parameter_square, sizeof(double), &parameter_matrix_bytes))
        return NUMELYX_INFERENCE_SIZE_OVERFLOW;

    weighted_moments = (double *) malloc(moment_bytes);
    if (gauss_newton_hessian != NULL)
        weighted_jacobian = (double *) malloc(jacobian_bytes);
    if (gradient != NULL)
        gradient_candidate = (double *) malloc(parameter_bytes);
    if (gauss_newton_hessian != NULL)
        hessian_candidate = (double *) malloc(parameter_matrix_bytes);
    if (weighted_moments == NULL ||
        (gauss_newton_hessian != NULL && weighted_jacobian == NULL) ||
        (gradient != NULL && gradient_candidate == NULL) ||
        (gauss_newton_hessian != NULL && hessian_candidate == NULL)) {
        status = NUMELYX_INFERENCE_OUT_OF_MEMORY;
        goto cleanup;
    }

    for (i = 0; i < moment_count; i++) {
        if (!isfinite(moments[i])) {
            status = NUMELYX_INFERENCE_NON_FINITE_INPUT;
            goto cleanup;
        }
        for (j = 0; j < moment_count; j++) {
            double value = weighting_matrix[i * moment_count + j];
            if (!isfinite(value)) {
                status = NUMELYX_INFERENCE_NON_FINITE_INPUT;
                goto cleanup;
            }
            if (fabs(value) > weighting_scale)
                weighting_scale = fabs(value);
            if (fabs(value - weighting_matrix[j * moment_count + i]) >
                1e-10 * fmax(1.0, weighting_scale)) {
                status = NUMELYX_INFERENCE_INVALID_ARGUMENT;
                goto cleanup;
            }
        }
        for (j = 0; j < parameter_count; j++) {
            if (!isfinite(moment_jacobian[i * parameter_count + j])) {
                status = NUMELYX_INFERENCE_NON_FINITE_INPUT;
                goto cleanup;
            }
        }
    }

    for (i = 0; i < moment_count; i++) {
        double value = 0.0;
        double correction = 0.0;
        for (j = 0; j < moment_count; j++) {
            gmm_kahan_add(&value, &correction,
                weighting_matrix[i * moment_count + j] * moments[j]);
        }
        if (!isfinite(value)) {
            status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
            goto cleanup;
        }
        weighted_moments[i] = value;
        gmm_kahan_add(&objective_value, &objective_correction,
            moments[i] * value);
        objective_absolute_sum += fabs(moments[i] * value);
    }
    if (!isfinite(objective_value) || !isfinite(objective_absolute_sum)) {
        status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
        goto cleanup;
    }
    if (objective_value < 0.0) {
        double tolerance = 64.0 * 2.2204460492503131e-16 *
            fmax(1.0, objective_absolute_sum);
        if (objective_value < -tolerance) {
            status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
            goto cleanup;
        }
        objective_value = 0.0;
    }

    if (gradient != NULL) {
        for (j = 0; j < parameter_count; j++) {
            double value = 0.0;
            double correction = 0.0;
            for (i = 0; i < moment_count; i++) {
                gmm_kahan_add(&value, &correction,
                    moment_jacobian[i * parameter_count + j] *
                    weighted_moments[i]);
            }
            gradient_candidate[j] = 2.0 * value;
            if (!isfinite(gradient_candidate[j])) {
                status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
                goto cleanup;
            }
        }
    }

    if (gauss_newton_hessian != NULL) {
        for (i = 0; i < moment_count; i++) {
            for (j = 0; j < parameter_count; j++) {
                double value = 0.0;
                double correction = 0.0;
                for (k = 0; k < moment_count; k++) {
                    gmm_kahan_add(&value, &correction,
                        weighting_matrix[i * moment_count + k] *
                        moment_jacobian[k * parameter_count + j]);
                }
                weighted_jacobian[i * parameter_count + j] = value;
            }
        }
        for (i = 0; i < parameter_count; i++) {
            for (j = 0; j < parameter_count; j++) {
                double value = 0.0;
                double correction = 0.0;
                for (k = 0; k < moment_count; k++) {
                    gmm_kahan_add(&value, &correction,
                        moment_jacobian[k * parameter_count + i] *
                        weighted_jacobian[k * parameter_count + j]);
                }
                hessian_candidate[i * parameter_count + j] = 2.0 * value;
                if (!isfinite(hessian_candidate[i * parameter_count + j])) {
                    status = NUMELYX_INFERENCE_NUMERICAL_FAILURE;
                    goto cleanup;
                }
            }
        }
        for (i = 0; i < parameter_count; i++) {
            for (j = i + 1; j < parameter_count; j++) {
                double symmetric =
                    0.5 * hessian_candidate[i * parameter_count + j] +
                    0.5 * hessian_candidate[j * parameter_count + i];
                hessian_candidate[i * parameter_count + j] = symmetric;
                hessian_candidate[j * parameter_count + i] = symmetric;
            }
        }
    }

    *objective = objective_value;
    if (gradient != NULL)
        memcpy(gradient, gradient_candidate, parameter_bytes);
    if (gauss_newton_hessian != NULL)
        memcpy(gauss_newton_hessian, hessian_candidate, parameter_matrix_bytes);
    status = NUMELYX_INFERENCE_SUCCESS;

cleanup:
    free(weighted_moments);
    free(weighted_jacobian);
    free(gradient_candidate);
    free(hessian_candidate);
    return status;
}


numelyx_core_inference_status numelyx_core_gmm_overidentification_log_p_value(
    double j_statistic,
    size_t moment_count,
    size_t parameter_count,
    double *log_p_value
)
{
    double value;
    numelyx_core_numeric_status status;
    if (log_p_value == NULL || moment_count == 0 || parameter_count == 0)
        return NUMELYX_INFERENCE_INVALID_ARGUMENT;
    if (!isfinite(j_statistic))
        return NUMELYX_INFERENCE_NON_FINITE_INPUT;
    if (j_statistic < 0.0 || moment_count <= parameter_count)
        return NUMELYX_INFERENCE_INVALID_ARGUMENT;
    status = numelyx_core_chi_square_log_survival(
        j_statistic, (double) (moment_count - parameter_count), &value);
    if (status != NUMERUS_NUMERIC_SUCCESS)
        return NUMELYX_INFERENCE_NUMERICAL_FAILURE;
    *log_p_value = value;
    return NUMELYX_INFERENCE_SUCCESS;
}
