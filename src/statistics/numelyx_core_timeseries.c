/**
 * @file numelyx_core_timeseries.c
 * @brief Stable autocovariance, autocorrelation, ACF, and differencing.
 */
#include "numelyx_core_timeseries.h"
#include "numelyx_core_size.h"

#include <math.h>
#include <stdlib.h>

static void timeseries_kahan_add(double *sum, double *correction, double value)
{
    double adjusted = value - *correction;
    double next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

static numelyx_core_timeseries_status timeseries_mean(
    const double *series, size_t count, double *mean
)
{
    size_t i;
    double sum = 0.0;
    double correction = 0.0;
    for (i = 0; i < count; i++) {
        if (!isfinite(series[i]))
            return NUMELYX_TIMESERIES_NON_FINITE_INPUT;
        timeseries_kahan_add(&sum, &correction, series[i]);
    }
    *mean = sum / (double) count;
    return isfinite(*mean)
        ? NUMELYX_TIMESERIES_SUCCESS
        : NUMELYX_TIMESERIES_NUMERICAL_FAILURE;
}

static numelyx_core_timeseries_status autocovariance_with_mean(
    const double *series, size_t count, size_t lag, double mean,
    bool unbiased, double *result
)
{
    size_t i;
    double sum = 0.0;
    double correction = 0.0;
    double denominator = (double) (unbiased ? count - lag : count);

    for (i = lag; i < count; i++) {
        double current = series[i] - mean;
        double lagged = series[i - lag] - mean;
        double term = current * lagged;
        if (!isfinite(current) || !isfinite(lagged) || !isfinite(term))
            return NUMELYX_TIMESERIES_NUMERICAL_FAILURE;
        timeseries_kahan_add(&sum, &correction, term);
    }
    *result = sum / denominator;
    return isfinite(*result)
        ? NUMELYX_TIMESERIES_SUCCESS
        : NUMELYX_TIMESERIES_NUMERICAL_FAILURE;
}

numelyx_core_timeseries_status numelyx_core_autocovariance(
    const double *series, size_t count, size_t lag,
    bool unbiased, double *result
)
{
    double mean;
    double value;
    numelyx_core_timeseries_status status;
    if (series == NULL || result == NULL || count == 0)
        return NUMELYX_TIMESERIES_INVALID_ARGUMENT;
    if (lag >= count) return NUMELYX_TIMESERIES_DOMAIN_ERROR;
    status = timeseries_mean(series, count, &mean);
    if (status != NUMELYX_TIMESERIES_SUCCESS) return status;
    status = autocovariance_with_mean(
        series, count, lag, mean, unbiased, &value);
    if (status != NUMELYX_TIMESERIES_SUCCESS) return status;
    *result = value;
    return NUMELYX_TIMESERIES_SUCCESS;
}

numelyx_core_timeseries_status numelyx_core_autocorrelation(
    const double *series, size_t count, size_t lag, double *result
)
{
    double mean;
    double variance;
    double covariance;
    numelyx_core_timeseries_status status;
    if (series == NULL || result == NULL || count == 0)
        return NUMELYX_TIMESERIES_INVALID_ARGUMENT;
    if (lag >= count) return NUMELYX_TIMESERIES_DOMAIN_ERROR;
    status = timeseries_mean(series, count, &mean);
    if (status != NUMELYX_TIMESERIES_SUCCESS) return status;
    status = autocovariance_with_mean(series, count, 0, mean, false, &variance);
    if (status != NUMELYX_TIMESERIES_SUCCESS) return status;
    if (variance <= 0.0) return NUMELYX_TIMESERIES_DOMAIN_ERROR;
    status = autocovariance_with_mean(
        series, count, lag, mean, false, &covariance);
    if (status != NUMELYX_TIMESERIES_SUCCESS) return status;
    *result = covariance / variance;
    return isfinite(*result)
        ? NUMELYX_TIMESERIES_SUCCESS
        : NUMELYX_TIMESERIES_NUMERICAL_FAILURE;
}

numelyx_core_timeseries_status numelyx_core_acf(
    const double *series, size_t count, size_t max_lag, double *acf_values
)
{
    size_t output_count;
    size_t output_bytes;
    size_t lag;
    double mean;
    double variance;
    double *candidate;
    numelyx_core_timeseries_status status;

    if (series == NULL || acf_values == NULL || count == 0)
        return NUMELYX_TIMESERIES_INVALID_ARGUMENT;
    if (max_lag >= count) return NUMELYX_TIMESERIES_DOMAIN_ERROR;
    if (max_lag == (size_t) -1)
        return NUMELYX_TIMESERIES_SIZE_OVERFLOW;
    output_count = max_lag + 1;
    if (!numelyx_core_size_multiply(
            output_count, sizeof(*candidate), &output_bytes))
        return NUMELYX_TIMESERIES_SIZE_OVERFLOW;
    status = timeseries_mean(series, count, &mean);
    if (status != NUMELYX_TIMESERIES_SUCCESS) return status;
    status = autocovariance_with_mean(series, count, 0, mean, false, &variance);
    if (status != NUMELYX_TIMESERIES_SUCCESS) return status;
    if (variance <= 0.0) return NUMELYX_TIMESERIES_DOMAIN_ERROR;

    candidate = (double *) malloc(output_bytes);
    if (candidate == NULL) return NUMELYX_TIMESERIES_OUT_OF_MEMORY;
    candidate[0] = 1.0;
    for (lag = 1; lag <= max_lag; lag++) {
        double covariance;
        status = autocovariance_with_mean(
            series, count, lag, mean, false, &covariance);
        if (status != NUMELYX_TIMESERIES_SUCCESS) {
            free(candidate);
            return status;
        }
        candidate[lag] = covariance / variance;
        if (!isfinite(candidate[lag])) {
            free(candidate);
            return NUMELYX_TIMESERIES_NUMERICAL_FAILURE;
        }
    }
    for (lag = 0; lag < output_count; lag++)
        acf_values[lag] = candidate[lag];
    free(candidate);
    return NUMELYX_TIMESERIES_SUCCESS;
}

numelyx_core_timeseries_status numelyx_core_difference(
    const double *series, size_t count, size_t lag,
    double *output, size_t *output_count
)
{
    size_t index;
    size_t result_count;
    if (series == NULL || output == NULL || count == 0)
        return NUMELYX_TIMESERIES_INVALID_ARGUMENT;
    if (lag == 0 || lag >= count)
        return NUMELYX_TIMESERIES_DOMAIN_ERROR;
    result_count = count - lag;
    for (index = 0; index < count; index++) {
        if (!isfinite(series[index]))
            return NUMELYX_TIMESERIES_NON_FINITE_INPUT;
    }
    for (index = 0; index < result_count; index++) {
        double difference = series[index + lag] - series[index];
        if (!isfinite(difference))
            return NUMELYX_TIMESERIES_NUMERICAL_FAILURE;
    }
    for (index = 0; index < result_count; index++)
        output[index] = series[index + lag] - series[index];
    if (output_count != NULL) *output_count = result_count;
    return NUMELYX_TIMESERIES_SUCCESS;
}
