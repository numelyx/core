/**
 * @file numelyx_core_special.c
 * @brief Regularized incomplete gamma and beta functions.
 *
 * Algorithms use a power series or modified Lentz continued fraction,
 * selecting the direct tail to avoid subtracting a small probability from 1.
 */
#include "numelyx_core_special.h"

#include <float.h>
#include <math.h>

#define NUMELYX_SPECIAL_MAX_ITERATIONS 100000
#define NUMELYX_SPECIAL_EPSILON (8.0 * DBL_EPSILON)
#define NUMELYX_SPECIAL_FPMIN (DBL_MIN / DBL_EPSILON)

static double clamp_probability(double value)
{
    if (value < 0.0 && value > -32.0 * DBL_EPSILON) return 0.0;
    if (value > 1.0 && value < 1.0 + 32.0 * DBL_EPSILON) return 1.0;
    return value;
}

/* Computes the continued fraction for the incomplete beta function. */
static bool beta_continued_fraction(
    double a, double b, double x, double *fraction
)
{
    size_t iteration;
    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    double h;

    if (fabs(d) < NUMELYX_SPECIAL_FPMIN) d = NUMELYX_SPECIAL_FPMIN;
    d = 1.0 / d;
    h = d;

    for (iteration = 1; iteration <= NUMELYX_SPECIAL_MAX_ITERATIONS; iteration++) {
        double m = (double) iteration;
        double m2 = 2.0 * m;
        double aa = m * (b - m) * x /
            ((qam + m2) * (a + m2));

        d = 1.0 + aa * d;
        if (fabs(d) < NUMELYX_SPECIAL_FPMIN) d = NUMELYX_SPECIAL_FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < NUMELYX_SPECIAL_FPMIN) c = NUMELYX_SPECIAL_FPMIN;
        d = 1.0 / d;
        h *= d * c;

        aa = -(a + m) * (qab + m) * x /
            ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < NUMELYX_SPECIAL_FPMIN) d = NUMELYX_SPECIAL_FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < NUMELYX_SPECIAL_FPMIN) c = NUMELYX_SPECIAL_FPMIN;
        d = 1.0 / d;
        {
            double delta = d * c;
            h *= delta;
            if (!isfinite(h) || !isfinite(delta)) return false;
            if (fabs(delta - 1.0) <= NUMELYX_SPECIAL_EPSILON) {
                *fraction = h;
                return true;
            }
        }
    }
    return false;
}

/* Lower regularized gamma series, convergent for x < a + 1. */
static bool gamma_p_series(double a, double x, double *value)
{
    size_t iteration;
    double sum = 1.0 / a;
    double term = sum;
    double ap = a;
    double log_scale = a * log(x) - x - lgamma(a);
    double scale = exp(log_scale);

    if (!isfinite(scale)) return false;
    for (iteration = 0; iteration < NUMELYX_SPECIAL_MAX_ITERATIONS; iteration++) {
        ap += 1.0;
        term *= x / ap;
        sum += term;
        if (!isfinite(sum) || !isfinite(term)) return false;
        if (fabs(term) <= fabs(sum) * NUMELYX_SPECIAL_EPSILON) {
            *value = clamp_probability(sum * scale);
            return isfinite(*value) && *value >= 0.0 && *value <= 1.0;
        }
    }
    return false;
}

/* Upper regularized gamma continued fraction, convergent for x >= a + 1. */
static bool gamma_q_fraction(double a, double x, double *value)
{
    size_t iteration;
    double b = x + 1.0 - a;
    double c = 1.0 / NUMELYX_SPECIAL_FPMIN;
    double d = 1.0 / (fabs(b) < NUMELYX_SPECIAL_FPMIN
        ? copysign(NUMELYX_SPECIAL_FPMIN, b == 0.0 ? 1.0 : b) : b);
    double h = d;
    double scale = exp(a * log(x) - x - lgamma(a));

    if (!isfinite(scale)) return false;
    for (iteration = 1; iteration <= NUMELYX_SPECIAL_MAX_ITERATIONS; iteration++) {
        double i = (double) iteration;
        double an = -i * (i - a);
        double delta;

        b += 2.0;
        d = an * d + b;
        if (fabs(d) < NUMELYX_SPECIAL_FPMIN) d = NUMELYX_SPECIAL_FPMIN;
        c = b + an / c;
        if (fabs(c) < NUMELYX_SPECIAL_FPMIN) c = NUMELYX_SPECIAL_FPMIN;
        d = 1.0 / d;
        delta = d * c;
        h *= delta;
        if (!isfinite(h) || !isfinite(delta)) return false;
        if (fabs(delta - 1.0) <= NUMELYX_SPECIAL_EPSILON) {
            *value = clamp_probability(scale * h);
            return isfinite(*value) && *value >= 0.0 && *value <= 1.0;
        }
    }
    return false;
}

static numelyx_core_numeric_status validate_gamma(
    double a, double x, double *result
)
{
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(a) || !isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    if (a <= 0.0 || x < 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_regularized_gamma_p(
    double a, double x, double *result
)
{
    double value;
    numelyx_core_numeric_status status = validate_gamma(a, x, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (x == 0.0) value = 0.0;
    else if (x < a + 1.0) {
        if (!gamma_p_series(a, x, &value))
            return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    } else {
        double q;
        if (!gamma_q_fraction(a, x, &q))
            return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
        value = 1.0 - q;
    }
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_regularized_gamma_q(
    double a, double x, double *result
)
{
    double value;
    numelyx_core_numeric_status status = validate_gamma(a, x, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (x == 0.0) value = 1.0;
    else if (x < a + 1.0) {
        double p;
        if (!gamma_p_series(a, x, &p))
            return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
        value = 1.0 - p;
    } else {
        if (!gamma_q_fraction(a, x, &value))
            return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_regularized_beta(
    double a, double b, double x, double *result
)
{
    double log_front;
    double value;
    double fraction;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(a) || !isfinite(b) || !isfinite(x))
        return NUMERUS_NUMERIC_NON_FINITE;
    if (a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    if (x == 0.0) value = 0.0;
    else if (x == 1.0) value = 1.0;
    else {
        log_front = lgamma(a + b) - lgamma(a) - lgamma(b) +
            a * log(x) + b * log1p(-x);
        if (!isfinite(log_front))
            return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
        if (x < (a + 1.0) / (a + b + 2.0)) {
            double log_tail;
            if (!beta_continued_fraction(a, b, x, &fraction) ||
                fraction <= 0.0)
                return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
            log_tail = log_front + log(fraction) - log(a);
            value = exp(log_tail);
        } else {
            double log_complement;
            if (!beta_continued_fraction(b, a, 1.0 - x, &fraction) ||
                fraction <= 0.0)
                return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
            log_complement = log_front + log(fraction) - log(b);
            value = 1.0 - exp(log_complement);
        }
        value = clamp_probability(value);
        if (!isfinite(value) || value < 0.0 || value > 1.0)
            return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}


#define NUMELYX_SPECIAL_INVERSE_MAX_ITERATIONS 2048

static numelyx_core_numeric_status validate_probability_inverse(
    double probability, double *result
)
{
    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(probability)) return NUMERUS_NUMERIC_NON_FINITE;
    if (probability <= 0.0 || probability >= 1.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    return NUMERUS_NUMERIC_SUCCESS;
}

static double safe_midpoint(double lower, double upper)
{
    return lower + (upper - lower) * 0.5;
}

numelyx_core_numeric_status numelyx_core_inverse_regularized_gamma_p(
    double a, double probability, double *result
)
{
    double lower = 0.0;
    double upper;
    double value;
    size_t iteration;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(a)) return NUMERUS_NUMERIC_NON_FINITE;
    if (a <= 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;
    status = validate_probability_inverse(probability, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (probability > 0.5) {
        return numelyx_core_inverse_regularized_gamma_q(
            a, 1.0 - probability, result);
    }

    upper = a > 1.0 ? a : 1.0;
    for (;;) {
        status = numelyx_core_regularized_gamma_p(a, upper, &value);
        if (status != NUMERUS_NUMERIC_SUCCESS)
            return status;
        if (value >= probability) break;
        if (upper > DBL_MAX * 0.5) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
        upper *= 2.0;
    }

    for (iteration = 0; iteration < NUMELYX_SPECIAL_INVERSE_MAX_ITERATIONS;
         iteration++) {
        double midpoint = safe_midpoint(lower, upper);
        if (midpoint == lower || midpoint == upper) break;
        status = numelyx_core_regularized_gamma_p(a, midpoint, &value);
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
        if (value < probability) lower = midpoint;
        else upper = midpoint;
    }
    *result = safe_midpoint(lower, upper);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_inverse_regularized_gamma_q(
    double a, double probability, double *result
)
{
    double lower = 0.0;
    double upper;
    double value;
    size_t iteration;
    numelyx_core_numeric_status status;

    if (!isfinite(a)) return NUMERUS_NUMERIC_NON_FINITE;
    if (a <= 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;
    status = validate_probability_inverse(probability, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (probability > 0.5) {
        return numelyx_core_inverse_regularized_gamma_p(
            a, 1.0 - probability, result);
    }

    upper = a > 1.0 ? a : 1.0;
    for (;;) {
        status = numelyx_core_regularized_gamma_q(a, upper, &value);
        if (status != NUMERUS_NUMERIC_SUCCESS)
            return status;
        if (value <= probability) break;
        if (upper > DBL_MAX * 0.5) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
        upper *= 2.0;
    }

    for (iteration = 0; iteration < NUMELYX_SPECIAL_INVERSE_MAX_ITERATIONS;
         iteration++) {
        double midpoint = safe_midpoint(lower, upper);
        if (midpoint == lower || midpoint == upper) break;
        status = numelyx_core_regularized_gamma_q(a, midpoint, &value);
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
        if (value > probability) lower = midpoint;
        else upper = midpoint;
    }
    *result = safe_midpoint(lower, upper);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_inverse_regularized_beta(
    double a, double b, double probability, double *result
)
{
    double lower = 0.0;
    double upper = 1.0;
    double target = probability;
    bool reflected = false;
    double value;
    size_t iteration;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(a) || !isfinite(b)) return NUMERUS_NUMERIC_NON_FINITE;
    if (a <= 0.0 || b <= 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;
    status = validate_probability_inverse(probability, result);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;

    if (probability > 0.5) {
        double temporary = a;
        a = b;
        b = temporary;
        target = 1.0 - probability;
        reflected = true;
    }

    for (iteration = 0; iteration < NUMELYX_SPECIAL_INVERSE_MAX_ITERATIONS;
         iteration++) {
        double midpoint = safe_midpoint(lower, upper);
        if (midpoint == lower || midpoint == upper) break;
        status = numelyx_core_regularized_beta(a, b, midpoint, &value);
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
        if (value < target) lower = midpoint;
        else upper = midpoint;
    }

    value = safe_midpoint(lower, upper);
    *result = reflected ? 1.0 - value : value;
    return NUMERUS_NUMERIC_SUCCESS;
}

/*
 * Digamma and trigamma use recurrence to shift x into the asymptotic region,
 * followed by Bernoulli-number expansions. This avoids reflection near poles;
 * the current contract is restricted to positive real arguments.
 */
numelyx_core_numeric_status numelyx_core_digamma(double x, double *result)
{
    double value = 0.0;
    double inverse;
    double inverse2;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    if (x <= 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;

    while (x < 8.0) {
        value -= 1.0 / x;
        x += 1.0;
        if (!isfinite(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    inverse = 1.0 / x;
    inverse2 = inverse * inverse;
    value += log(x) - 0.5 * inverse
        - inverse2 * (1.0 / 12.0
        - inverse2 * (1.0 / 120.0
        - inverse2 * (1.0 / 252.0
        - inverse2 * (1.0 / 240.0
        - inverse2 * (5.0 / 660.0
        - inverse2 * (691.0 / 32760.0
        - inverse2 * (1.0 / 12.0)))))));

    if (!isfinite(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_trigamma(double x, double *result)
{
    double value = 0.0;
    double inverse;
    double inverse2;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    if (x <= 0.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;

    while (x < 8.0) {
        double term = 1.0 / x;
        term *= term;
        value += term;
        x += 1.0;
        if (!isfinite(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    inverse = 1.0 / x;
    inverse2 = inverse * inverse;
    value += inverse + 0.5 * inverse2
        + inverse * inverse2 * (1.0 / 6.0
        - inverse2 * (1.0 / 30.0
        - inverse2 * (1.0 / 42.0
        - inverse2 * (1.0 / 30.0
        - inverse2 * (5.0 / 66.0
        - inverse2 * (691.0 / 2730.0
        - inverse2 * (7.0 / 6.0)))))));

    if (!isfinite(value)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}


numelyx_core_numeric_status numelyx_core_erf(double x, double *result)
{
    double gamma_value;
    double square;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    if (x == 0.0) {
        *result = x;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (fabs(x) > sqrt(DBL_MAX)) {
        *result = copysign(1.0, x);
        return NUMERUS_NUMERIC_SUCCESS;
    }

    square = x * x;
    status = numelyx_core_regularized_gamma_p(0.5, square, &gamma_value);
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    *result = copysign(gamma_value, x);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_erfc(double x, double *result)
{
    double gamma_value;
    double square;
    double value;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(x)) return NUMERUS_NUMERIC_NON_FINITE;
    if (x == 0.0) {
        *result = 1.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (fabs(x) > sqrt(DBL_MAX)) {
        *result = x > 0.0 ? 0.0 : 2.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }

    square = x * x;
    if (x > 0.0) {
        status = numelyx_core_regularized_gamma_q(0.5, square, &gamma_value);
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
        value = gamma_value;
    } else {
        status = numelyx_core_regularized_gamma_p(0.5, square, &gamma_value);
        if (status != NUMERUS_NUMERIC_SUCCESS) return status;
        value = 1.0 + gamma_value;
    }
    if (!isfinite(value) || value < 0.0 || value > 2.0)
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}


numelyx_core_numeric_status numelyx_core_inverse_erf(double y, double *result)
{
    double magnitude;
    double gamma_argument;
    double root;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(y)) return NUMERUS_NUMERIC_NON_FINITE;
    if (y < -1.0 || y > 1.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;
    if (y == 0.0) {
        *result = y;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (y == 1.0 || y == -1.0) {
        *result = copysign(INFINITY, y);
        return NUMERUS_NUMERIC_SUCCESS;
    }

    magnitude = fabs(y);
    /* The first term of the inverse-erf series is correctly rounded here. */
    if (magnitude < 1e-8) {
        *result = copysign(magnitude * sqrt(acos(-1.0)) * 0.5, y);
        return NUMERUS_NUMERIC_SUCCESS;
    }

    if (magnitude <= 0.5) {
        status = numelyx_core_inverse_regularized_gamma_p(
            0.5, magnitude, &gamma_argument);
    } else {
        status = numelyx_core_inverse_regularized_gamma_q(
            0.5, 1.0 - magnitude, &gamma_argument);
    }
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    root = sqrt(gamma_argument);
    if (!isfinite(root)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = copysign(root, y);
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_inverse_erfc(double y, double *result)
{
    double gamma_argument;
    double magnitude;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(y)) return NUMERUS_NUMERIC_NON_FINITE;
    if (y < 0.0 || y > 2.0) return NUMERUS_NUMERIC_DOMAIN_ERROR;
    if (y == 0.0) {
        *result = INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (y == 2.0) {
        *result = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (y == 1.0) {
        *result = 0.0;
        return NUMERUS_NUMERIC_SUCCESS;
    }

    if (y < 1.0) {
        status = numelyx_core_inverse_regularized_gamma_q(
            0.5, y, &gamma_argument);
        magnitude = 1.0;
    } else {
        status = numelyx_core_inverse_regularized_gamma_p(
            0.5, 2.0 - y, &gamma_argument);
        magnitude = -1.0;
    }
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    gamma_argument = sqrt(gamma_argument);
    if (!isfinite(gamma_argument))
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = magnitude * gamma_argument;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_normal_quantile(
    double probability, double *result
)
{
    double tail_probability;
    double quantile;
    numelyx_core_numeric_status status;

    if (result == NULL) return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    if (!isfinite(probability)) return NUMERUS_NUMERIC_NON_FINITE;
    if (probability < 0.0 || probability > 1.0)
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    if (probability == 0.0) {
        *result = -INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }
    if (probability == 1.0) {
        *result = INFINITY;
        return NUMERUS_NUMERIC_SUCCESS;
    }

    if (probability <= 0.5) {
        tail_probability = 2.0 * probability;
        status = numelyx_core_inverse_erfc(tail_probability, &quantile);
        quantile = -sqrt(2.0) * quantile;
    } else {
        tail_probability = 2.0 * (1.0 - probability);
        status = numelyx_core_inverse_erfc(tail_probability, &quantile);
        quantile *= sqrt(2.0);
    }
    if (status != NUMERUS_NUMERIC_SUCCESS) return status;
    if (!isfinite(quantile)) return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    *result = quantile;
    return NUMERUS_NUMERIC_SUCCESS;
}
