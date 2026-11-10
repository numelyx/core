/**
 * @file numelyx_core_fuzzy.h
 * @brief Composable scalar fuzzy-membership and fuzzy-operator primitives.
 */
#ifndef NUMELYX_CORE_FUZZY_H
#define NUMELYX_CORE_FUZZY_H

typedef enum {
    NUMELYX_FUZZY_SUCCESS = 0,
    NUMELYX_FUZZY_INVALID_ARGUMENT,
    NUMELYX_FUZZY_NON_FINITE_INPUT,
    NUMELYX_FUZZY_DOMAIN_ERROR
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

#endif
