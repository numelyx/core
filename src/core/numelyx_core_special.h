/**
 * @file numelyx_core_special.h
 * @brief Stable special-function primitives used by probability algorithms.
 */
#ifndef NUMELYX_CORE_SPECIAL_H
#define NUMELYX_CORE_SPECIAL_H

#include "numelyx_core_numeric.h"

/**
 * Regularized lower incomplete gamma P(a,x), for finite a > 0 and x >= 0.
 * The result is in [0,1] and is written only on success.
 */
numelyx_core_numeric_status numelyx_core_regularized_gamma_p(
    double a, double x, double *result
);

/**
 * Regularized upper incomplete gamma Q(a,x), evaluated through a direct
 * continued fraction in its numerically favorable region.
 */
numelyx_core_numeric_status numelyx_core_regularized_gamma_q(
    double a, double x, double *result
);

/**
 * Regularized incomplete beta I_x(a,b), for finite a,b > 0 and x in [0,1].
 */
numelyx_core_numeric_status numelyx_core_regularized_beta(
    double a, double b, double x, double *result
);

#endif
