/**
 * @file numelyx_core_numeric.c
 * @brief Scale-aware floating-point comparison and tolerance helpers.
 *
 * @details This implementation follows the focused subsystem headers and
 * preserves their ownership, validation, and status-reporting contracts.
 */

#include "numelyx_core_numeric.h"

#include <math.h>

numelyx_core_numeric_status numelyx_core_numeric_sigmoid(
    double input,
    double *result
)
{
    double value;

    if (result == NULL) {
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    }
    if (!isfinite(input)) {
        return NUMERUS_NUMERIC_NON_FINITE;
    }

    if (input >= 0.0) {
        value = 1.0 / (1.0 + exp(-input));
    } else {
        double exponential = exp(input);
        value = exponential / (1.0 + exponential);
    }

    if (!isfinite(value)) {
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_numeric_logit(
    double probability,
    double *result
)
{
    double value;

    if (result == NULL) {
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    }
    if (!isfinite(probability)) {
        return NUMERUS_NUMERIC_NON_FINITE;
    }
    if (probability <= 0.0 || probability >= 1.0) {
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    }

    value = log(probability) - log1p(-probability);
    if (!isfinite(value)) {
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_numeric_softplus(
    double input,
    double *result
)
{
    double value;

    if (result == NULL) {
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    }
    if (!isfinite(input)) {
        return NUMERUS_NUMERIC_NON_FINITE;
    }

    value = fmax(input, 0.0) + log1p(exp(-fabs(input)));
    if (!isfinite(value)) {
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_numeric_log_sigmoid(
    double input,
    double *result
)
{
    double value;

    if (result == NULL) {
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    }
    if (!isfinite(input)) {
        return NUMERUS_NUMERIC_NON_FINITE;
    }

    if (input >= 0.0) {
        value = -log1p(exp(-input));
    } else {
        value = input - log1p(exp(input));
    }

    if (!isfinite(value)) {
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_numeric_log_sum_exp(
    const double *values,
    size_t count,
    double *result
)
{
    size_t index;
    double maximum;
    double sum = 0.0;
    double compensation = 0.0;
    double value;

    if (values == NULL || result == NULL || count == 0) {
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    }

    maximum = values[0];
    if (!isfinite(maximum)) {
        return NUMERUS_NUMERIC_NON_FINITE;
    }

    for (index = 1; index < count; index++) {
        if (!isfinite(values[index])) {
            return NUMERUS_NUMERIC_NON_FINITE;
        }
        if (values[index] > maximum) {
            maximum = values[index];
        }
    }

    for (index = 0; index < count; index++) {
        double term = exp(values[index] - maximum);
        double adjusted;
        double next_sum;

        if (!isfinite(term)) {
            return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
        }

        /* Kahan summation is inexpensive here and limits accumulation error. */
        adjusted = term - compensation;
        next_sum = sum + adjusted;
        if (!isfinite(next_sum)) {
            return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
        }
        compensation = (next_sum - sum) - adjusted;
        sum = next_sum;
    }

    value = maximum + log(sum);
    if (!isfinite(value)) {
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}


/*
 * Stable log(1 - exp(x)). The split at -log(2) avoids catastrophic
 * cancellation both near zero and when exp(x) is small.
 */
numelyx_core_numeric_status numelyx_core_numeric_log1mexp(
    double input,
    double *result
)
{
    const double negative_log_two = -0.693147180559945309417232121458176568;
    double value;

    if (result == NULL) {
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    }
    if (!isfinite(input)) {
        return NUMERUS_NUMERIC_NON_FINITE;
    }
    if (input > 0.0) {
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    }

    if (input == 0.0) {
        value = -INFINITY;
    } else if (input < negative_log_two) {
        value = log1p(-exp(input));
    } else {
        value = log(-expm1(input));
    }

    /* -INFINITY is the exact and documented value at the endpoint. */
    if (isnan(value) || value == INFINITY) {
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}

numelyx_core_numeric_status numelyx_core_numeric_log_diff_exp(
    double a,
    double b,
    double *result
)
{
    double difference;
    double correction;
    double value;
    numelyx_core_numeric_status status;

    if (result == NULL) {
        return NUMERUS_NUMERIC_INVALID_ARGUMENT;
    }
    if (!isfinite(a) || !isfinite(b)) {
        return NUMERUS_NUMERIC_NON_FINITE;
    }
    if (b > a) {
        return NUMERUS_NUMERIC_DOMAIN_ERROR;
    }

    difference = b - a;
    if (difference == 0.0) {
        value = -INFINITY;
    } else {
        status = numelyx_core_numeric_log1mexp(difference, &correction);
        if (status != NUMERUS_NUMERIC_SUCCESS) {
            return status;
        }
        value = a + correction;
    }

    /* -INFINITY is valid only when the two logarithms are equal. */
    if (isnan(value) || value == INFINITY ||
        (value == -INFINITY && a != b)) {
        return NUMERUS_NUMERIC_NUMERICAL_FAILURE;
    }

    *result = value;
    return NUMERUS_NUMERIC_SUCCESS;
}
