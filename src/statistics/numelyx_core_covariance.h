/**
 * @file numelyx_core_covariance.h
 * @brief Generic sandwich covariance estimators for score/estimating equations.
 */
#ifndef NUMELYX_CORE_COVARIANCE_H
#define NUMELYX_CORE_COVARIANCE_H

#include <stddef.h>

typedef enum {
    NUMELYX_COVARIANCE_SUCCESS = 0,
    NUMELYX_COVARIANCE_INVALID_ARGUMENT,
    NUMELYX_COVARIANCE_NON_FINITE_INPUT,
    NUMELYX_COVARIANCE_DOMAIN_ERROR,
    NUMELYX_COVARIANCE_SIZE_OVERFLOW,
    NUMELYX_COVARIANCE_OUT_OF_MEMORY,
    NUMELYX_COVARIANCE_NUMERICAL_FAILURE
} numelyx_core_covariance_status;

typedef enum {
    NUMELYX_COVARIANCE_HC0 = 0,
    NUMELYX_COVARIANCE_HC1
} numelyx_core_covariance_correction;

typedef enum {
    NUMELYX_COVARIANCE_CLUSTER_CR0 = 0,
    NUMELYX_COVARIANCE_CLUSTER_CR1
} numelyx_core_cluster_correction;

/**
 * Compute B (sum_i s_i s_i^T) B^T from row-major observation scores.
 * HC1 applies n/(n-p), requiring n>p. Output is unchanged on failure.
 */
numelyx_core_covariance_status numelyx_core_sandwich_covariance(
    const double *bread,
    const double *scores,
    size_t observation_count,
    size_t parameter_count,
    numelyx_core_covariance_correction correction,
    double *covariance
);

/**
 * Cluster-robust sandwich using dense cluster IDs in [0, cluster_count).
 * CR1 applies G/(G-1); CR0 is uncorrected.
 */
numelyx_core_covariance_status numelyx_core_cluster_sandwich_covariance(
    const double *bread,
    const double *scores,
    const size_t *cluster_ids,
    size_t observation_count,
    size_t parameter_count,
    size_t cluster_count,
    numelyx_core_cluster_correction correction,
    double *covariance
);

/**
 * Newey-West HAC sandwich with Bartlett weights and the requested lag
 * bandwidth (clamped to observation_count-1).
 */
numelyx_core_covariance_status numelyx_core_hac_sandwich_covariance(
    const double *bread,
    const double *scores,
    size_t observation_count,
    size_t parameter_count,
    size_t bandwidth,
    double *covariance
);

typedef enum {
    NUMELYX_COVARIANCE_HC_TYPE_0 = 0,
    NUMELYX_COVARIANCE_HC_TYPE_1,
    NUMELYX_COVARIANCE_HC_TYPE_2,
    NUMELYX_COVARIANCE_HC_TYPE_3
} numelyx_core_hc_type;

/**
 * Regression-score HC0-HC3 covariance from X, residuals, and bread.
 * Leverage is required for HC2/HC3 and may be NULL for HC0/HC1.
 */
numelyx_core_covariance_status numelyx_core_hc_sandwich_covariance(
    const double *bread,
    const double *design,
    const double *residuals,
    const double *leverage,
    size_t observation_count,
    size_t parameter_count,
    numelyx_core_hc_type type,
    double *covariance
);

#endif
