/**
 * @file numelyx_core_fuzzy.c
 * @brief Allocation-free fuzzy membership and algebraic operator primitives.
 */
#include "numelyx_core_fuzzy.h"

#include <float.h>
#include <stddef.h>
#include <math.h>
#include <stdint.h>

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


static numelyx_core_fuzzy_status validate_intervals(
    double left_lower, double left_upper,
    double right_lower, double right_upper,
    double *result_lower, double *result_upper
)
{
    if (result_lower == NULL || result_upper == NULL ||
        result_lower == result_upper)
        return NUMELYX_FUZZY_INVALID_ARGUMENT;
    if (!isfinite(left_lower) || !isfinite(left_upper) ||
        !isfinite(right_lower) || !isfinite(right_upper))
        return NUMELYX_FUZZY_NON_FINITE_INPUT;
    if (left_lower > left_upper || right_lower > right_upper)
        return NUMELYX_FUZZY_DOMAIN_ERROR;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_interval_add(
    double left_lower, double left_upper, double right_lower, double right_upper,
    double *result_lower, double *result_upper
)
{
    double lower;
    double upper;
    numelyx_core_fuzzy_status status = validate_intervals(
        left_lower, left_upper, right_lower, right_upper,
        result_lower, result_upper);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    lower = left_lower + right_lower;
    upper = left_upper + right_upper;
    if (!isfinite(lower) || !isfinite(upper))
        return NUMELYX_FUZZY_NUMERICAL_FAILURE;
    *result_lower = lower;
    *result_upper = upper;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_interval_subtract(
    double left_lower, double left_upper, double right_lower, double right_upper,
    double *result_lower, double *result_upper
)
{
    double lower;
    double upper;
    numelyx_core_fuzzy_status status = validate_intervals(
        left_lower, left_upper, right_lower, right_upper,
        result_lower, result_upper);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    lower = left_lower - right_upper;
    upper = left_upper - right_lower;
    if (!isfinite(lower) || !isfinite(upper))
        return NUMELYX_FUZZY_NUMERICAL_FAILURE;
    *result_lower = lower;
    *result_upper = upper;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_interval_multiply(
    double left_lower, double left_upper, double right_lower, double right_upper,
    double *result_lower, double *result_upper
)
{
    double products[4];
    double lower;
    double upper;
    size_t index;
    numelyx_core_fuzzy_status status = validate_intervals(
        left_lower, left_upper, right_lower, right_upper,
        result_lower, result_upper);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;

    products[0] = left_lower * right_lower;
    products[1] = left_lower * right_upper;
    products[2] = left_upper * right_lower;
    products[3] = left_upper * right_upper;
    lower = products[0];
    upper = products[0];
    for (index = 0; index < 4; index++) {
        if (!isfinite(products[index]))
            return NUMELYX_FUZZY_NUMERICAL_FAILURE;
        if (products[index] < lower) lower = products[index];
        if (products[index] > upper) upper = products[index];
    }
    *result_lower = lower;
    *result_upper = upper;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_interval_divide(
    double numerator_lower, double numerator_upper,
    double denominator_lower, double denominator_upper,
    double *result_lower, double *result_upper
)
{
    double reciprocal_lower;
    double reciprocal_upper;
    double lower;
    double upper;
    numelyx_core_fuzzy_status status = validate_intervals(
        numerator_lower, numerator_upper,
        denominator_lower, denominator_upper,
        result_lower, result_upper);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;
    if (denominator_lower <= 0.0 && denominator_upper >= 0.0)
        return NUMELYX_FUZZY_DOMAIN_ERROR;

    reciprocal_lower = 1.0 / denominator_upper;
    reciprocal_upper = 1.0 / denominator_lower;
    if (reciprocal_lower > reciprocal_upper) {
        double temporary = reciprocal_lower;
        reciprocal_lower = reciprocal_upper;
        reciprocal_upper = temporary;
    }

    {
        double products[4] = {
            numerator_lower * reciprocal_lower,
            numerator_lower * reciprocal_upper,
            numerator_upper * reciprocal_lower,
            numerator_upper * reciprocal_upper
        };
        size_t index;
        lower = products[0];
        upper = products[0];
        for (index = 0; index < 4; index++) {
            if (!isfinite(products[index]))
                return NUMELYX_FUZZY_NUMERICAL_FAILURE;
            if (products[index] < lower) lower = products[index];
            if (products[index] > upper) upper = products[index];
        }
    }
    *result_lower = lower;
    *result_upper = upper;
    return NUMELYX_FUZZY_SUCCESS;
}


static numelyx_core_fuzzy_status validate_sampled_fuzzy_set(
    const double *x_values, const double *memberships,
    size_t count, double *result, double *coordinate_scale
)
{
    size_t index;
    double scale = 0.0;
    if (x_values == NULL || memberships == NULL || result == NULL || count == 0)
        return NUMELYX_FUZZY_INVALID_ARGUMENT;
    for (index = 0; index < count; index++) {
        if (!isfinite(x_values[index]) || !isfinite(memberships[index]))
            return NUMELYX_FUZZY_NON_FINITE_INPUT;
        if (memberships[index] < 0.0 || memberships[index] > 1.0)
            return NUMELYX_FUZZY_DOMAIN_ERROR;
        if (fabs(x_values[index]) > scale) scale = fabs(x_values[index]);
    }
    *coordinate_scale = scale == 0.0 ? 1.0 : scale;
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_defuzzify_centroid(
    const double *x_values, const double *memberships, size_t count,
    double *result
)
{
    size_t index;
    double coordinate_scale;
    double total_weight = 0.0;
    double weighted_sum = 0.0;
    double correction = 0.0;
    numelyx_core_fuzzy_status status = validate_sampled_fuzzy_set(
        x_values, memberships, count, result, &coordinate_scale);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;

    for (index = 0; index < count; index++) {
        double term;
        double adjusted;
        double next;
        total_weight += memberships[index];
        term = (x_values[index] / coordinate_scale) * memberships[index];
        adjusted = term - correction;
        next = weighted_sum + adjusted;
        correction = (next - weighted_sum) - adjusted;
        weighted_sum = next;
    }
    if (total_weight <= 0.0)
        return NUMELYX_FUZZY_DOMAIN_ERROR;
    {
        double value = (weighted_sum / total_weight) * coordinate_scale;
        if (!isfinite(value)) return NUMELYX_FUZZY_NUMERICAL_FAILURE;
        *result = value;
    }
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_defuzzify_mean_of_maxima(
    const double *x_values, const double *memberships, size_t count,
    double *result
)
{
    size_t index;
    size_t maxima_count = 0;
    double maximum = 0.0;
    double coordinate_scale;
    double scaled_sum = 0.0;
    double correction = 0.0;
    numelyx_core_fuzzy_status status = validate_sampled_fuzzy_set(
        x_values, memberships, count, result, &coordinate_scale);
    if (status != NUMELYX_FUZZY_SUCCESS) return status;

    for (index = 0; index < count; index++) {
        if (memberships[index] > maximum) maximum = memberships[index];
    }
    if (maximum <= 0.0) return NUMELYX_FUZZY_DOMAIN_ERROR;

    for (index = 0; index < count; index++) {
        if (memberships[index] == maximum) {
            double term = x_values[index] / coordinate_scale;
            double adjusted = term - correction;
            double next = scaled_sum + adjusted;
            correction = (next - scaled_sum) - adjusted;
            scaled_sum = next;
            maxima_count++;
        }
    }
    {
        double value = (scaled_sum / (double) maxima_count) * coordinate_scale;
        if (!isfinite(value)) return NUMELYX_FUZZY_NUMERICAL_FAILURE;
        *result = value;
    }
    return NUMELYX_FUZZY_SUCCESS;
}


numelyx_core_fuzzy_status numelyx_core_fuzzy_mamdani_aggregate(
    const double *rule_strengths,
    const double *consequent_memberships,
    size_t rule_count,
    size_t sample_count,
    numelyx_core_fuzzy_implication implication,
    numelyx_core_fuzzy_aggregation aggregation,
    double *output_memberships
)
{
    size_t matrix_elements;
    size_t rule;
    size_t sample;

    if (rule_strengths == NULL || consequent_memberships == NULL ||
        output_memberships == NULL || rule_count == 0 || sample_count == 0)
        return NUMELYX_FUZZY_INVALID_ARGUMENT;
    if ((implication != NUMELYX_FUZZY_IMPLICATION_MIN &&
         implication != NUMELYX_FUZZY_IMPLICATION_PRODUCT) ||
        (aggregation != NUMELYX_FUZZY_AGGREGATION_MAX &&
         aggregation != NUMELYX_FUZZY_AGGREGATION_PROBABILISTIC_SUM &&
         aggregation != NUMELYX_FUZZY_AGGREGATION_BOUNDED_SUM))
        return NUMELYX_FUZZY_INVALID_ARGUMENT;
    if (rule_count > SIZE_MAX / sample_count)
        return NUMELYX_FUZZY_NUMERICAL_FAILURE;
    matrix_elements = rule_count * sample_count;

    for (rule = 0; rule < rule_count; rule++) {
        if (!isfinite(rule_strengths[rule]))
            return NUMELYX_FUZZY_NON_FINITE_INPUT;
        if (rule_strengths[rule] < 0.0 || rule_strengths[rule] > 1.0)
            return NUMELYX_FUZZY_DOMAIN_ERROR;
    }
    for (sample = 0; sample < matrix_elements; sample++) {
        if (!isfinite(consequent_memberships[sample]))
            return NUMELYX_FUZZY_NON_FINITE_INPUT;
        if (consequent_memberships[sample] < 0.0 ||
            consequent_memberships[sample] > 1.0)
            return NUMELYX_FUZZY_DOMAIN_ERROR;
    }

    for (sample = 0; sample < sample_count; sample++) {
        double aggregate = 0.0;
        for (rule = 0; rule < rule_count; rule++) {
            double consequent =
                consequent_memberships[rule * sample_count + sample];
            double implied = implication == NUMELYX_FUZZY_IMPLICATION_MIN
                ? fmin(rule_strengths[rule], consequent)
                : rule_strengths[rule] * consequent;

            if (aggregation == NUMELYX_FUZZY_AGGREGATION_MAX) {
                if (implied > aggregate) aggregate = implied;
            } else if (aggregation ==
                NUMELYX_FUZZY_AGGREGATION_PROBABILISTIC_SUM) {
                aggregate = aggregate + implied - aggregate * implied;
            } else {
                aggregate = fmin(1.0, aggregate + implied);
            }
        }
        output_memberships[sample] = aggregate;
    }
    return NUMELYX_FUZZY_SUCCESS;
}

numelyx_core_fuzzy_status numelyx_core_fuzzy_sugeno_aggregate(
    const double *rule_strengths,
    const double *rule_outputs,
    size_t rule_count,
    double *output
)
{
    size_t index;
    double total_strength = 0.0;
    double output_scale = 0.0;
    double weighted_sum = 0.0;
    double correction = 0.0;

    if (rule_strengths == NULL || rule_outputs == NULL ||
        output == NULL || rule_count == 0)
        return NUMELYX_FUZZY_INVALID_ARGUMENT;

    for (index = 0; index < rule_count; index++) {
        if (!isfinite(rule_strengths[index]) || !isfinite(rule_outputs[index]))
            return NUMELYX_FUZZY_NON_FINITE_INPUT;
        if (rule_strengths[index] < 0.0 || rule_strengths[index] > 1.0)
            return NUMELYX_FUZZY_DOMAIN_ERROR;
        total_strength += rule_strengths[index];
        if (fabs(rule_outputs[index]) > output_scale)
            output_scale = fabs(rule_outputs[index]);
    }
    if (total_strength <= 0.0) return NUMELYX_FUZZY_DOMAIN_ERROR;
    if (output_scale == 0.0) {
        *output = 0.0;
        return NUMELYX_FUZZY_SUCCESS;
    }

    for (index = 0; index < rule_count; index++) {
        double term = (rule_outputs[index] / output_scale) *
            rule_strengths[index];
        double adjusted = term - correction;
        double next = weighted_sum + adjusted;
        correction = (next - weighted_sum) - adjusted;
        weighted_sum = next;
    }
    {
        double value = (weighted_sum / total_strength) * output_scale;
        if (!isfinite(value)) return NUMELYX_FUZZY_NUMERICAL_FAILURE;
        *output = value;
    }
    return NUMELYX_FUZZY_SUCCESS;
}
