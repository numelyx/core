/**
 * @file numelyx_core_fuzzy.c
 * @brief Allocation-free fuzzy membership and algebraic operator primitives.
 */
#include "numelyx_core_fuzzy.h"

#include <float.h>
#include <math.h>

static numelyx_core_fuzzy_status validate_membership_input(
    double x, double *membership
)
{
    if (membership == NULL) return NUMELYX_FUZZY_INVALID_ARGUMENT;
    if (!isfinite(x)) return NUMELYX_FUZZY_NON_FINITE_INPUT;
    return NUMELYX_FUZZY_SUCCESS;
}

static numelyx_core_fuzzy_status validate_unit_memberships(
    double left, double right, double *result
)
{
    if (result == NULL) return NUMELYX_FUZZY_INVALID_ARGUMENT;
    if (!isfinite(left) || !isfinite(right))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (left < 0.0 || left > 1.0 || right < 0.0 || right > 1.0)
        return NUMELYX_FUZZY_DOMAIN_ERROR;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_triangular_membership(
    double x, double a, double b, double c, double *membership
)
{
    double value;
    numelyx_core_fuzzy_status status = validate_membership_input(x, membership);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    if (!isfinite(a) || !isfinite(b) || !isfinite(c))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (!(a <= b && b <= c) || a == c)
        return NUMELYX_FUZZY_DOMAIN_ERROR;

    if ((a == b && x == a) || (b == c && x == c)) value = 1.0;
    else if (x <= a || x >= c) value = 0.0;
    else if (x == b) value = 1.0;
    else if (x < b) value = (x - a) / (b - a);
    else value = (c - x) / (c - b);

    *membership = value;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_trapezoidal_membership(
    double x, double a, double b, double c, double d, double *membership
)
{
    double value;
    numelyx_core_fuzzy_status status = validate_membership_input(x, membership);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    if (!isfinite(a) || !isfinite(b) || !isfinite(c) || !isfinite(d))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (!(a <= b && b <= c && c <= d) || a == d)
        return NUMELYX_FUZZY_DOMAIN_ERROR;

    if ((a == b && x == a) || (c == d && x == d)) value = 1.0;
    else if (x <= a || x >= d) value = 0.0;
    else if (x >= b && x <= c) value = 1.0;
    else if (x < b) value = (x - a) / (b - a);
    else value = (d - x) / (d - c);

    *membership = value;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_gaussian_membership(
    double x, double center, double sigma, double *membership
)
{
    double difference;
    double z;
    double value;
    numelyx_core_fuzzy_status status = validate_membership_input(x, membership);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    if (!isfinite(center) || !isfinite(sigma))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (sigma <= 0.0) return NUMELYX_FUZZY_DOMAIN_ERROR;

    difference = x - center;
    if (!isfinite(difference)) {
        value = 0.0;
    } else {
        z = difference / sigma;
        if (!isfinite(z) || fabs(z) > sqrt(DBL_MAX) * sqrt(2.0)) {
            value = 0.0;
        } else {
            value = exp(-0.5 * z * z);
        }
    }
    *membership = value;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_bell_membership(
    double x, double width, double shape, double center, double *membership
)
{
    double difference;
    double log_term;
    double value;
    numelyx_core_fuzzy_status status = validate_membership_input(x, membership);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    if (!isfinite(width) || !isfinite(shape) || !isfinite(center))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (width <= 0.0 || shape <= 0.0)
        return NUMELYX_FUZZY_DOMAIN_ERROR;

    difference = x - center;
    if (difference == 0.0) {
        value = 1.0;
    } else if (!isfinite(difference)) {
        value = 0.0;
    } else {
        log_term = 2.0 * shape * (log(fabs(difference)) - log(width));
        if (log_term > 0.0) {
            double inverse = exp(-log_term);
            value = inverse / (1.0 + inverse);
        } else {
            value = 1.0 / (1.0 + exp(log_term));
        }
    }
    *membership = value;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_sigmoid_membership(
    double x, double slope, double center, double *membership
)
{
    double difference;
    double argument;
    double value;
    numelyx_core_fuzzy_status status = validate_membership_input(x, membership);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    if (!isfinite(slope) || !isfinite(center))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (slope == 0.0) {
        *membership = 0.5;
        return NUMELYX_FUZZY_SUCCESS;
    }

    difference = x - center;
    if (!isfinite(difference)) {
        argument = ((x > center) == (slope > 0.0)) ? INFINITY : -INFINITY;
    } else {
        argument = slope * difference;
    }

    if (argument >= 0.0) {
        value = 1.0 / (1.0 + exp(-argument));
    } else {
        double exponential = exp(argument);
        value = exponential / (1.0 + exponential);
    }
    *membership = value;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_complement(
    double membership, double *result
)
{
    numelyx_core_fuzzy_status status =
        validate_unit_memberships(membership, 0.0, result);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    *result = 1.0 - membership;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_tnorm_min(
    double left, double right, double *result
)
{
    numelyx_core_fuzzy_status status =
        validate_unit_memberships(left, right, result);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    *result = left < right ? left : right;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_tnorm_product(
    double left, double right, double *result
)
{
    numelyx_core_fuzzy_status status =
        validate_unit_memberships(left, right, result);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    *result = left * right;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_tnorm_lukasiewicz(
    double left, double right, double *result
)
{
    numelyx_core_fuzzy_status status =
        validate_unit_memberships(left, right, result);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    *result = fmax(0.0, left + right - 1.0);
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_snorm_max(
    double left, double right, double *result
)
{
    numelyx_core_fuzzy_status status =
        validate_unit_memberships(left, right, result);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    *result = left > right ? left : right;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_snorm_probabilistic_sum(
    double left, double right, double *result
)
{
    numelyx_core_fuzzy_status status =
        validate_unit_memberships(left, right, result);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    *result = left + right - left * right;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_snorm_bounded_sum(
    double left, double right, double *result
)
{
    numelyx_core_fuzzy_status status =
        validate_unit_memberships(left, right, result);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    *result = fmin(1.0, left + right);
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_triangular_alpha_cut(
    double a, double b, double c, double alpha,
    double *lower, double *upper
)
{
    if (lower == NULL || upper == NULL)
        return NUMELYX_FUZZY_INVALID_ARGUMENT;
    if (!isfinite(a) || !isfinite(b) || !isfinite(c) || !isfinite(alpha))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (!(a <= b && b <= c) || a == c || alpha < 0.0 || alpha > 1.0)
        return NUMELYX_FUZZY_DOMAIN_ERROR;
    *lower = (1.0 - alpha) * a + alpha * b;
    *upper = (1.0 - alpha) * c + alpha * b;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_trapezoidal_alpha_cut(
    double a, double b, double c, double d, double alpha,
    double *lower, double *upper
)
{
    if (lower == NULL || upper == NULL)
        return NUMELYX_FUZZY_INVALID_ARGUMENT;
    if (!isfinite(a) || !isfinite(b) || !isfinite(c) ||
        !isfinite(d) || !isfinite(alpha))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (!(a <= b && b <= c && c <= d) || a == d ||
        alpha < 0.0 || alpha > 1.0)
        return NUMELYX_FUZZY_DOMAIN_ERROR;
    *lower = (1.0 - alpha) * a + alpha * b;
    *upper = (1.0 - alpha) * d + alpha * c;
    return NUMELYX_FUZZY_SUCCESS;
}
