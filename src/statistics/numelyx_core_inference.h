/**
 * @file numelyx_core_inference.h
 * @brief Generic delta-method and Wald inference primitives.
 */
#ifndef NUMELYX_CORE_INFERENCE_H
#define NUMELYX_CORE_INFERENCE_H

#include <stddef.h>

typedef enum {
    NUMELYX_INFERENCE_SUCCESS = 0,
    NUMELYX_INFERENCE_INVALID_ARGUMENT,
    NUMELYX_INFERENCE_NON_FINITE_INPUT,
    NUMELYX_INFERENCE_SIZE_OVERFLOW,
    NUMELYX_INFERENCE_OUT_OF_MEMORY,
    NUMELYX_INFERENCE_SINGULAR_COVARIANCE,
    NUMELYX_INFERENCE_NUMERICAL_FAILURE
} numelyx_core_inference_status;

/** Compute g^T V g from a row-major covariance matrix V. */
numelyx_core_inference_status numelyx_core_delta_method_variance(
    const double *gradient,
    const double *covariance,
    size_t parameter_count,
    double *variance
);

/** Compute the delta-method standard error sqrt(g^T V g). */
numelyx_core_inference_status numelyx_core_delta_method_standard_error(
    const double *gradient,
    const double *covariance,
    size_t parameter_count,
    double *standard_error
);

/**
 * Wald test for R*estimate = null_value.
 *
 * Computes W = d^T (R V R^T)^(-1) d using Cholesky solves, where
 * d = R*estimate-null_value. log_p_value is the chi-square upper-tail log
 * probability with restriction_count degrees of freedom.
 */
numelyx_core_inference_status numelyx_core_wald_test(
    const double *estimate,
    const double *null_value,
    const double *restriction_matrix,
    const double *covariance,
    size_t parameter_count,
    size_t restriction_count,
    double *wald_statistic,
    double *log_p_value
);

#endif
