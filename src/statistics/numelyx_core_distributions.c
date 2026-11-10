/**
 * @file numelyx_core_distributions.c
 * @brief Stable scalar normal-distribution primitives.
 */
#include "numelyx_core_distributions.h"
#include "numelyx_core_special.h"
#include "numelyx_core_rng.h"

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
    double log_q;
    double gamma_argument;
    double log_small_tail;
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

    if (fabs(z) > sqrt(DBL_MAX) * sqrt(2.0) * sqrt(2.0)) {
        log_q = -INFINITY;
    } else {
        gamma_argument = 0.5 * z * z;
        status = numelyx_core_regularized_gamma_log_q(
            0.5, gamma_argument, &log_q
        );
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    }

    log_small_tail = NUMELYX_LOG_HALF + log_q;
    if (z < 0.0) {
        *log_cdf = log_small_tail;
        *log_survival = log1p(-exp(log_small_tail));
    } else {
        *log_cdf = log1p(-exp(log_small_tail));
        *log_survival = log_small_tail;
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


static numelyx_core_numeric_status validate_student_t(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
)
{
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(x) || !isfinite(degrees_of_freedom) ||
        !isfinite(location) || !isfinite(scale))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (degrees_of_freedom <= 0.0 || scale <= 0.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    return NUMERUS_NUMERIC_SUCCESS;
}

static double log1pexp_stable(double value)
{
    if (value > 0.0) return value + log1p(exp(-value));
    return log1p(exp(value));
}

numelyx_core_numeric_status numelyx_core_student_t_log_pdf(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
)
{
    const double log_pi = 1.1447298858494001741434273513530587116;
    double z;
    double log_kernel;
    double value;
    numelyx_core_numeric_status status =
        validate_student_t(x, degrees_of_freedom, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;

    z = standardized_difference(x, location, scale);
    if (!isfinite(z)) {
        *result = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (z == 0.0) {
        log_kernel = 0.0;
    } else {
        double log_ratio = 2.0 * log(fabs(z)) - log(degrees_of_freedom);
        log_kernel = log1pexp_stable(log_ratio);
    }
    value = lgamma(0.5 * (degrees_of_freedom + 1.0)) -
        lgamma(0.5 * degrees_of_freedom) -
        0.5 * (log(degrees_of_freedom) + log_pi) -
        log(scale) -
        0.5 * (degrees_of_freedom + 1.0) * log_kernel;
    if (isnan(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_student_t_pdf(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
)
{
    double log_density;
    numelyx_core_numeric_status status;
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    status = numelyx_core_student_t_log_pdf(
        x, degrees_of_freedom, location, scale, &log_density);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_density);
    return NUMERUS_NUMERIC_SUCCESS;
}

static numelyx_core_numeric_status student_t_log_tails(
    double x, double degrees_of_freedom, double location, double scale,
    double *log_cdf, double *log_survival
)
{
    double z = standardized_difference(x, location, scale);
    double log_tail;
    double ratio;
    double beta_x;
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

    if (z == 0.0) {
        log_tail = NUMELYX_LOG_HALF;
    } else {
        if (fabs(z) > sqrt(degrees_of_freedom)) {
            ratio = sqrt(degrees_of_freedom) / fabs(z);
            beta_x = (ratio * ratio) / (1.0 + ratio * ratio);
        } else {
            ratio = fabs(z) / sqrt(degrees_of_freedom);
            beta_x = 1.0 / (1.0 + ratio * ratio);
        }
        status = numelyx_core_regularized_beta_log_lower(
            0.5 * degrees_of_freedom, 0.5, beta_x, &log_tail);
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
        log_tail += NUMELYX_LOG_HALF;
    }

    if (z < 0.0) {
        *log_cdf = log_tail;
        *log_survival = log1p(-exp(log_tail));
    } else {
        *log_cdf = log1p(-exp(log_tail));
        *log_survival = log_tail;
    }
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_student_t_log_cdf(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_student_t(x, degrees_of_freedom, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = student_t_log_tails(
        x, degrees_of_freedom, location, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_cdf;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_student_t_log_survival(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_student_t(x, degrees_of_freedom, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = student_t_log_tails(
        x, degrees_of_freedom, location, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_survival;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_student_t_cdf(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_student_t(x, degrees_of_freedom, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = student_t_log_tails(
        x, degrees_of_freedom, location, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_cdf);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_student_t_survival(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_student_t(x, degrees_of_freedom, location, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = student_t_log_tails(
        x, degrees_of_freedom, location, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_survival);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_student_t_quantile(
    double probability, double degrees_of_freedom, double location,
    double scale, double *result
)
{
    double tail_probability;
    double beta_x;
    double log_magnitude;
    double magnitude;
    double value;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(probability) || !isfinite(degrees_of_freedom) ||
        !isfinite(location) || !isfinite(scale))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (probability < 0.0 || probability > 1.0 ||
        degrees_of_freedom <= 0.0 || scale <= 0.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    if (probability == 0.0) {
        *result = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (probability == 1.0) {
        *result = INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (probability == 0.5) {
        *result = location;
        return NUMERUS_NUMERIC_SUCCESS;
    }

    if (probability < 0.5) {
        tail_probability = 2.0 * probability;
    } else {
        tail_probability = 2.0 * (1.0 - probability);
    }
    status = numelyx_core_inverse_regularized_beta(
        0.5 * degrees_of_freedom, 0.5, tail_probability, &beta_x);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (beta_x <= 0.0) {
        magnitude = INFINITY;
    } else {
        log_magnitude = 0.5 * (
            log(degrees_of_freedom) + log1p(-beta_x) - log(beta_x));
        magnitude = exp(log_magnitude);
    }
    value = location + (probability < 0.5 ? -1.0 : 1.0) * scale * magnitude;
    if (isnan(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

static numelyx_core_numeric_status validate_chi_square(
    double degrees_of_freedom, double *result
)
{
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(degrees_of_freedom)) return NUMERUS_NUMERIC_NON_FINITE;
    if (degrees_of_freedom <= 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_chi_square_log_pdf(
    double x, double degrees_of_freedom, double *result
)
{
    double shape;
    double value;
    numelyx_core_numeric_status status =
        validate_chi_square(degrees_of_freedom, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    if (x < 0.0) {
        *result = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    shape = 0.5 * degrees_of_freedom;
    if (x == 0.0) {
        if (shape < 1.0) value = INFINITY;
        else if (shape == 1.0) value = -log(2.0);
        else value = -INFINITY;
    } else {
        value = (shape - 1.0) * log(x) - 0.5 * x -
            shape * log(2.0) - lgamma(shape);
    }
    if (isnan(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_chi_square_pdf(
    double x, double degrees_of_freedom, double *result
)
{
    double log_density;
    numelyx_core_numeric_status status;
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    status = numelyx_core_chi_square_log_pdf(x, degrees_of_freedom, &log_density);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_density);
    return NUMERUS_NUMERIC_SUCCESS;
}

static numelyx_core_numeric_status chi_square_log_tails(
    double x, double degrees_of_freedom,
    double *log_cdf, double *log_survival
)
{
    double log_q;
    numelyx_core_numeric_status status;
    if (x < 0.0) {
        *log_cdf = -INFINITY;
        *log_survival = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (x == 0.0) {
        *log_cdf = -INFINITY;
        *log_survival = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    status = numelyx_core_regularized_gamma_log_q(
        0.5 * degrees_of_freedom, 0.5 * x, &log_q);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *log_survival = log_q;
    return numelyx_core_numeric_log1mexp(log_q, log_cdf);
}

numelyx_core_numeric_status numelyx_core_chi_square_log_cdf(
    double x, double degrees_of_freedom, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_chi_square(degrees_of_freedom, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    status = chi_square_log_tails(x, degrees_of_freedom, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_cdf;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_chi_square_log_survival(
    double x, double degrees_of_freedom, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_chi_square(degrees_of_freedom, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    status = chi_square_log_tails(x, degrees_of_freedom, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_survival;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_chi_square_cdf(
    double x, double degrees_of_freedom, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_chi_square(degrees_of_freedom, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    status = chi_square_log_tails(x, degrees_of_freedom, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_cdf);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_chi_square_survival(
    double x, double degrees_of_freedom, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_chi_square(degrees_of_freedom, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    status = chi_square_log_tails(x, degrees_of_freedom, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_survival);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_chi_square_quantile(
    double probability, double degrees_of_freedom, double *result
)
{
    double gamma_quantile;
    numelyx_core_numeric_status status;
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(probability) || !isfinite(degrees_of_freedom))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (probability < 0.0 || probability > 1.0 || degrees_of_freedom <= 0.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    if (probability == 0.0) {
        *result = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (probability == 1.0) {
        *result = INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    status = numelyx_core_inverse_regularized_gamma_p(
        0.5 * degrees_of_freedom, probability, &gamma_quantile);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = 2.0 * gamma_quantile;
    return NUMERUS_NUMERIC_SUCCESS;
}


numelyx_core_numeric_status numelyx_core_normal_sample(
    numelyx_core_rng *rng, double location, double scale, double *result
)
{
    numelyx_core_rng_status rng_status;

    if (rng == NULL || result == NULL)
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(location) || !isfinite(scale))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (scale <= 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;

    rng_status = numelyx_core_rng_normal_affine(rng, location, scale, result);
    if (rng_status == NUMERUS_RNG_SUCCESS) return NUMERUS_NUMERIC_SUCCESS;
    if (rng_status == NUMERUS_RNG_INVALID_ARGUMENT)
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
}


static numelyx_core_numeric_status validate_gamma_distribution(
    double x, double shape, double scale, double *result
)
{
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(x) || !isfinite(shape) || !isfinite(scale))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (shape <= 0.0 || scale <= 0.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_gamma_log_pdf(
    double x, double shape, double scale, double *result
)
{
    double value;
    numelyx_core_numeric_status status =
        validate_gamma_distribution(x, shape, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (x < 0.0) {
        *result = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (x == 0.0) {
        if (shape < 1.0) value = INFINITY;
        else if (shape == 1.0) value = -log(scale);
        else value = -INFINITY;
    } else {
        double scaled_x = x / scale;
        value = (shape - 1.0) * log(x) - scaled_x -
            lgamma(shape) - shape * log(scale);
    }
    if (isnan(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_gamma_pdf(
    double x, double shape, double scale, double *result
)
{
    double log_density;
    numelyx_core_numeric_status status;
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    status = numelyx_core_gamma_log_pdf(x, shape, scale, &log_density);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_density);
    return NUMERUS_NUMERIC_SUCCESS;
}

static numelyx_core_numeric_status gamma_distribution_log_tails(
    double x, double shape, double scale,
    double *log_cdf, double *log_survival
)
{
    double log_q;
    double scaled_x;
    numelyx_core_numeric_status status;

    if (x < 0.0) {
        *log_cdf = -INFINITY;
        *log_survival = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (x == 0.0) {
        *log_cdf = -INFINITY;
        *log_survival = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    scaled_x = x / scale;
    if (!isfinite(scaled_x)) {
        *log_cdf = 0.0;
        *log_survival = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    status = numelyx_core_regularized_gamma_log_q(shape, scaled_x, &log_q);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *log_survival = log_q;
    if (log_q == -INFINITY) {
        *log_cdf = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    return numelyx_core_numeric_log1mexp(log_q, log_cdf);
}

numelyx_core_numeric_status numelyx_core_gamma_log_cdf(
    double x, double shape, double scale, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_gamma_distribution(x, shape, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = gamma_distribution_log_tails(
        x, shape, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_cdf;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_gamma_log_survival(
    double x, double shape, double scale, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_gamma_distribution(x, shape, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = gamma_distribution_log_tails(
        x, shape, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_survival;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_gamma_cdf(
    double x, double shape, double scale, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_gamma_distribution(x, shape, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = gamma_distribution_log_tails(
        x, shape, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_cdf);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_gamma_survival(
    double x, double shape, double scale, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_gamma_distribution(x, shape, scale, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = gamma_distribution_log_tails(
        x, shape, scale, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_survival);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_gamma_quantile(
    double probability, double shape, double scale, double *result
)
{
    double gamma_quantile;
    numelyx_core_numeric_status status;
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(probability) || !isfinite(shape) || !isfinite(scale))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (probability < 0.0 || probability > 1.0 ||
        shape <= 0.0 || scale <= 0.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    if (probability == 0.0) {
        *result = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (probability == 1.0) {
        *result = INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    status = numelyx_core_inverse_regularized_gamma_p(
        shape, probability, &gamma_quantile);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = scale * gamma_quantile;
    return NUMERUS_NUMERIC_SUCCESS;
}

static numelyx_core_numeric_status validate_beta_distribution(
    double x, double alpha, double beta, double *result
)
{
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(x) || !isfinite(alpha) || !isfinite(beta))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (alpha <= 0.0 || beta <= 0.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_beta_log_pdf(
    double x, double alpha, double beta, double *result
)
{
    double value;
    double log_beta;
    numelyx_core_numeric_status status =
        validate_beta_distribution(x, alpha, beta, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (x < 0.0 || x > 1.0) {
        *result = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    log_beta = lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
    if (x == 0.0) {
        if (alpha < 1.0) value = INFINITY;
        else if (alpha == 1.0) value = -log_beta;
        else value = -INFINITY;
    } else if (x == 1.0) {
        if (beta < 1.0) value = INFINITY;
        else if (beta == 1.0) value = -log_beta;
        else value = -INFINITY;
    } else {
        value = (alpha - 1.0) * log(x) +
            (beta - 1.0) * log1p(-x) - log_beta;
    }
    if (isnan(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_beta_pdf(
    double x, double alpha, double beta, double *result
)
{
    double log_density;
    numelyx_core_numeric_status status;
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    status = numelyx_core_beta_log_pdf(x, alpha, beta, &log_density);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_density);
    return NUMERUS_NUMERIC_SUCCESS;
}

static numelyx_core_numeric_status beta_distribution_log_tails(
    double x, double alpha, double beta,
    double *log_cdf, double *log_survival
)
{
    if (x <= 0.0) {
        *log_cdf = -INFINITY;
        *log_survival = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (x >= 1.0) {
        *log_cdf = 0.0;
        *log_survival = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    {
        numelyx_core_numeric_status status =
            numelyx_core_regularized_beta_log_lower(
                alpha, beta, x, log_cdf);
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
        return numelyx_core_regularized_beta_log_upper(
            alpha, beta, x, log_survival);
    }
}

numelyx_core_numeric_status numelyx_core_beta_log_cdf(
    double x, double alpha, double beta, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_beta_distribution(x, alpha, beta, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = beta_distribution_log_tails(
        x, alpha, beta, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_cdf;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_beta_log_survival(
    double x, double alpha, double beta, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_beta_distribution(x, alpha, beta, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = beta_distribution_log_tails(
        x, alpha, beta, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = log_survival;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_beta_cdf(
    double x, double alpha, double beta, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_beta_distribution(x, alpha, beta, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = beta_distribution_log_tails(
        x, alpha, beta, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_cdf);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_beta_survival(
    double x, double alpha, double beta, double *result
)
{
    double log_cdf;
    double log_survival;
    numelyx_core_numeric_status status =
        validate_beta_distribution(x, alpha, beta, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    status = beta_distribution_log_tails(
        x, alpha, beta, &log_cdf, &log_survival);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = exp(log_survival);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_beta_quantile(
    double probability, double alpha, double beta, double *result
)
{
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(probability) || !isfinite(alpha) || !isfinite(beta))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (probability < 0.0 || probability > 1.0 ||
        alpha <= 0.0 || beta <= 0.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    if (probability == 0.0) {
        *result = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (probability == 1.0) {
        *result = 1.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    return numelyx_core_inverse_regularized_beta(
        alpha, beta, probability, result);
}
