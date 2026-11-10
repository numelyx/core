/**
 * @file numelyx_core_fuzzy.h
 * @brief Composable scalar fuzzy-membership and fuzzy-operator primitives.
 */
#ifndef NUMELYX_CORE_FUZZY_H
#define NUMELYX_CORE_FUZZY_H

#include <stddef.h>

typedef enum {
    NUMELYX_FUZZY_SUCCESS = 0,
    NUMELYX_FUZZY_INVALID_ARGUMENT,
    NUMELYX_FUZZY_NON_FINITE_INPUT,
    NUMELYX_FUZZY_DOMAIN_ERROR,
    NUMELYX_FUZZY_NUMERICAL_FAILURE
} numelyx_core_fuzzy_status;

numelyx_core_fuzzy_status numelyx_core_fuzzy_triangular_membership(
    double x, double a, double b, double c, double *membership
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_trapezoidal_membership(
    double x, double a, double b, double c, double d, double *membership
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_gaussian_membership(
    double x, double center, double sigma, double *membership
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_bell_membership(
    double x, double width, double shape, double center, double *membership
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_sigmoid_membership(
    double x, double slope, double center, double *membership
);

numelyx_core_fuzzy_status numelyx_core_fuzzy_complement(
    double membership, double *result
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_tnorm_min(
    double left, double right, double *result
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_tnorm_product(
    double left, double right, double *result
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_tnorm_lukasiewicz(
    double left, double right, double *result
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_snorm_max(
    double left, double right, double *result
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_snorm_probabilistic_sum(
    double left, double right, double *result
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_snorm_bounded_sum(
    double left, double right, double *result
);

numelyx_core_fuzzy_status numelyx_core_fuzzy_triangular_alpha_cut(
    double a, double b, double c, double alpha,
    double *lower, double *upper
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_trapezoidal_alpha_cut(
    double a, double b, double c, double d, double alpha,
    double *lower, double *upper
);

/** Interval arithmetic used by fuzzy alpha-cut extension operations. */
numelyx_core_fuzzy_status numelyx_core_fuzzy_interval_add(
    double left_lower, double left_upper, double right_lower, double right_upper,
    double *result_lower, double *result_upper
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_interval_subtract(
    double left_lower, double left_upper, double right_lower, double right_upper,
    double *result_lower, double *result_upper
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_interval_multiply(
    double left_lower, double left_upper, double right_lower, double right_upper,
    double *result_lower, double *result_upper
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_interval_divide(
    double numerator_lower, double numerator_upper,
    double denominator_lower, double denominator_upper,
    double *result_lower, double *result_upper
);

/** Discrete centroid and mean-of-maxima defuzzification for sampled sets. */
numelyx_core_fuzzy_status numelyx_core_fuzzy_defuzzify_centroid(
    const double *x_values, const double *memberships, size_t count,
    double *result
);
numelyx_core_fuzzy_status numelyx_core_fuzzy_defuzzify_mean_of_maxima(
    const double *x_values, const double *memberships, size_t count,
    double *result
);

typedef enum {
    NUMELYX_FUZZY_IMPLICATION_MIN = 0,
    NUMELYX_FUZZY_IMPLICATION_PRODUCT
} numelyx_core_fuzzy_implication;

typedef enum {
    NUMELYX_FUZZY_AGGREGATION_MAX = 0,
    NUMELYX_FUZZY_AGGREGATION_PROBABILISTIC_SUM,
    NUMELYX_FUZZY_AGGREGATION_BOUNDED_SUM
} numelyx_core_fuzzy_aggregation;

/**
 * Aggregate a rule-by-sample consequent-membership matrix using Mamdani
 * implication and the selected s-norm. Input and output buffers must not
 * overlap. Output is unchanged if validation fails.
 */
numelyx_core_fuzzy_status numelyx_core_fuzzy_mamdani_aggregate(
    const double *rule_strengths,
    const double *consequent_memberships,
    size_t rule_count,
    size_t sample_count,
    numelyx_core_fuzzy_implication implication,
    numelyx_core_fuzzy_aggregation aggregation,
    double *output_memberships
);

/** Weighted-average Sugeno aggregation; output is unchanged on failure. */
numelyx_core_fuzzy_status numelyx_core_fuzzy_sugeno_aggregate(
    const double *rule_strengths,
    const double *rule_outputs,
    size_t rule_count,
    double *output
);

#endif
