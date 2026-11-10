/**
 * @file numelyx_core_gmm.h
 * @brief Generic GMM moment-objective evaluation primitives.
 */
#ifndef NUMELYX_CORE_GMM_H
#define NUMELYX_CORE_GMM_H

#include <stddef.h>
#include "numelyx_core_inference.h"

/**
 * Evaluate Q = gᵀWg, gradient 2JᵀWg, and optional Gauss-Newton curvature
 * 2JᵀWJ. W is a symmetric q-by-q weighting matrix and J is a row-major
 * q-by-p moment Jacobian. Gradient and curvature outputs are optional.
 */
numelyx_core_inference_status numelyx_core_gmm_evaluate(
    const double *moments,
    const double *weighting_matrix,
    const double *moment_jacobian,
    size_t moment_count,
    size_t parameter_count,
    double *objective,
    double *gradient,
    double *gauss_newton_hessian
);

#endif
