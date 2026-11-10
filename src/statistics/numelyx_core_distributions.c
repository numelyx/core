/**
 * @file numelyx_core_distributions.c
 * @brief Stable scalar normal-distribution primitives.
 */
#include "numelyx_core_distributions.h"
#include "numelyx_core_special.h"

#include <float.h>
#include <math.h>

#define NUMELYX_LOG_SQRT_TWO_PI 0.9189385332046727417803297364056176399
#define NUMELYX_LOG_HALF  (-0.693147180559945309417232121458176568)

static numelyx_core_numeric_status validate_normal_arguments(
    double x, double location, double scale, double *result
)
{
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(x) || !isfinite(location) || !isfinite(scale))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (scale <= 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;
    return NUMERUS_NUMERIC_SUCCESS;
}

static double standardized_difference(double x, double location, double scale)
{
    double difference = x - location;
    if (!isfinite(difference)) {
        return x > location ? INFINITY : -INFINITY;
    }
    return difference / scale;
}

static numelyx_core_numeric_status normal_log_tails(
    double x,
    double location,
    double scale,
    double *log_cdf,
    double *log_survival
)
{
    double z = standardized_difference(x, location, scale);
    double q;
    double log_q;
    double gamma_argument;
    numelyx_core_numeric_status status;

    if (isinf(z)) {
        if (z < 0.0) {
            *log_cdf = -INFINITY;
            *log_survival = 0.0;
        } else {
            *log_cdf = 0.0;
            *log_survival = -INFINITY;
        }
        return NUMERUS_NUMERIC_SUCCESS;
    }

    if (fabs(z) > sqrt(DBL_MAX) * sqrt(2.0)) {
        q = 0.0;
    } else {
        gamma_argument = 0.5 * z * z;
        status = numelyx_core_regularized_gamma_q(
            0.5, gamma_argument, &q
        );
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    }

    log_q = q == 0.0 ? -INFINITY : log(q);
    if (z < 0.0) {
        *log_cdf = NUMELYX_LOG_HALF + log_q;
        *log_survival = log1p(-0.5 * q);
    } else {
        *log_cdf = log1p(-0.5 * q);
        *log_survival = NUMELYX_LOG_HALF + log_q;
    }
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_normal_log_pdf(
    double x, double location, double scale, double *result
)
{
    double z;
    double value;
    numelyx_core_numeric_status status =
        validate_normal_arguments(x, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;

    z = standardized_difference(x, location, scale);
    if (!isfinite(z) || fabs(z) > sqrt(DBL_MAX)) {
        *result = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    value = -0.5 * z * z - log(scale) - NUMELYX_LOG_SQRT_TWO_PI;
    if (isnan(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_normal_pdf(
    double x, double location, double scale, double *result
)
{
    double log_density;
    numelyx_core_numeric_status status;
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    status = numelyx_core_normal_log_pdf(x, location, scale, &log_density);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_density);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_normal_log_cdf(
    double x, double location, double scale, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_normal_arguments(x, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = normal_log_tails(x, location, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_cdf;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_normal_log_survival(
    double x, double location, double scale, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_normal_arguments(x, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = normal_log_tails(x, location, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_survival;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_normal_cdf(
    double x, double location, double scale, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_normal_arguments(x, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = normal_log_tails(x, location, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_cdf);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_normal_survival(
    double x, double location, double scale, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_normal_arguments(x, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = normal_log_tails(x, location, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_survival);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_normal_distribution_quantile(
    double probability, double location, double scale, double *result
)
{
    double standard_quantile;
    double value;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(probability) || !isfinite(location) || !isfinite(scale))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (probability < 0.0 || probability > 1.0 || scale <= 0.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;

    status = numelyx_core_normal_quantile(probability, &standard_quantile);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    value = location + scale * standard_quantile;
    if (isnan(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}
