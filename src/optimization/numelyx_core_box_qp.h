/**
 * @file numelyx_core_box_qp.h
 * @brief Internal solver for convex box-constrained quadratic programs.
 */
#ifndef NUMELYX_CORE_BOX_QP_H
#define NUMELYX_CORE_BOX_QP_H
#include <stddef.h>
typedef enum {
    NUMELYX_QP_CONVERGED = 0,
    NUMELYX_QP_ITERATION_LIMIT,
    NUMELYX_QP_UNBOUNDED,
    NUMELYX_QP_INVALID_ARGUMENT,
    NUMELYX_QP_NON_FINITE_INPUT,
    NUMELYX_QP_NOT_POSITIVE_SEMIDEFINITE,
    NUMELYX_QP_SIZE_OVERFLOW,
    NUMELYX_QP_OUT_OF_MEMORY,
    NUMELYX_QP_NUMERICAL_FAILURE
} numelyx_core_qp_status;
typedef struct {
    size_t max_iterations;
    double projected_gradient_tolerance;
    double positive_semidefinite_tolerance;
} numelyx_core_box_qp_options;
numelyx_core_qp_status numelyx_core_box_qp_default_options(
    numelyx_core_box_qp_options *options
);
numelyx_core_qp_status numelyx_core_box_qp_solve(
    const double *hessian, const double *linear_term,
    const double *lower_bounds, const double *upper_bounds,
    const double *initial_point, size_t variable_count,
    const numelyx_core_box_qp_options *options, double *solution,
    double *objective_value, size_t *iterations
);
#endif
