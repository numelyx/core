/**
 * @file numelyx_core_special.h
 * @brief Stable special-function primitives used by probability algorithms.
 */
#ifndef NUMELYX_CORE_SPECIAL_H
#define NUMELYX_CORE_SPECIAL_H

#include "numelyx_core_numeric.h"

/** Regularized lower incomplete gamma P(a,x), finite a > 0 and x >= 0. */
numelyx_core_numeric_status numelyx_core_regularized_gamma_p(
    double a, double x, double *result
);

/** Regularized upper incomplete gamma Q(a,x), finite a > 0 and x >= 0. */
numelyx_core_numeric_status numelyx_core_regularized_gamma_q(
    double a, double x, double *result
);

/** Regularized incomplete beta I_x(a,b), finite a,b > 0 and x in [0,1]. */
numelyx_core_numeric_status numelyx_core_regularized_beta(
    double a, double b, double x, double *result
);

/**
 * Inverse of P(a,x). Requires finite a > 0 and probability in (0,1).
 * Uses a bracketed monotone solve; output is unchanged on failure.
 */
numelyx_core_numeric_status numelyx_core_inverse_regularized_gamma_p(
    double a, double probability, double *result
);

/**
 * Inverse of Q(a,x). Requires finite a > 0 and probability in (0,1).
 * Evaluates the upper tail directly during inversion.
 */
numelyx_core_numeric_status numelyx_core_inverse_regularized_gamma_q(
    double a, double probability, double *result
);

/**
 * Inverse of I_x(a,b). Requires finite a,b > 0 and probability in (0,1).
 * Uses beta symmetry for upper-tail probabilities.
 */
numelyx_core_numeric_status numelyx_core_inverse_regularized_beta(
    double a, double b, double probability, double *result
);

/** Digamma psi(x), for finite x > 0. */
numelyx_core_numeric_status numelyx_core_digamma(double x, double *result);

/** Trigamma psi_1(x), for finite x > 0. */
numelyx_core_numeric_status numelyx_core_trigamma(double x, double *result);

#endif
