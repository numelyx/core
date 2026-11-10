/**
 * @file numelyx_core_timeseries.h
 * @brief Generic univariate time-series primitives.
 */
#ifndef NUMELYX_CORE_TIMESERIES_H
#define NUMELYX_CORE_TIMESERIES_H

#include <stdbool.h>
#include <stddef.h>

typedef enum {
    NUMELYX_TIMESERIES_SUCCESS = 0,
    NUMELYX_TIMESERIES_INVALID_ARGUMENT,
    NUMELYX_TIMESERIES_NON_FINITE_INPUT,
    NUMELYX_TIMESERIES_DOMAIN_ERROR,
    NUMELYX_TIMESERIES_SIZE_OVERFLOW,
    NUMELYX_TIMESERIES_OUT_OF_MEMORY,
    NUMELYX_TIMESERIES_NUMERICAL_FAILURE
} numelyx_core_timeseries_status;

/** Lag-k autocovariance about the full-sample mean; denominator n or n-k. */
numelyx_core_timeseries_status numelyx_core_autocovariance(
    const double *series, size_t observation_count, size_t lag,
    bool unbiased, double *result
);

/** Biased autocorrelation estimate rho(k)=gamma(k)/gamma(0). */
numelyx_core_timeseries_status numelyx_core_autocorrelation(
    const double *series, size_t observation_count, size_t lag,
    double *result
);

/** Compute ACF values for lags 0 through max_lag into a caller buffer. */
numelyx_core_timeseries_status numelyx_core_acf(
    const double *series, size_t observation_count, size_t max_lag,
    double *acf_values
);

/**
 * Apply one lag-difference: output[i] = series[i+lag] - series[i].
 * The output contains observation_count-lag values; output_count is optional.
 */
numelyx_core_timeseries_status numelyx_core_difference(
    const double *series, size_t observation_count, size_t lag,
    double *output, size_t *output_count
);

#endif
