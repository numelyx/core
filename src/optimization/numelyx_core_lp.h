/**
 * @file numelyx_core_lp.h
 * @brief Internal linear-programming solver for continuous LPs.
 */
#ifndef NUMELYX_CORE_LP_H
#define NUMELYX_CORE_LP_H

#include <stddef.h>

typedef enum {
    NUMELYX_LP_SUCCESS = 0,
    NUMELYX_LP_OPTIMAL,
    NUMELYX_LP_INFEASIBLE,
    NUMELYX_LP_UNBOUNDED,
    NUMELYX_LP_INVALID_ARGUMENT,
    NUMELYX_LP_NON_FINITE_INPUT,
    NUMELYX_LP_SIZE_OVERFLOW,
    NUMELYX_LP_OUT_OF_MEMORY,
    NUMELYX_LP_NUMERICAL_FAILURE,
    NUMELYX_LP_ITERATION_LIMIT
} numelyx_core_lp_status;

typedef struct {
    size_t max_iterations;
    double feasibility_tolerance;
    double optimality_tolerance;
} numelyx_core_lp_options;

/**
 * Solve max c^T x subject to A x <= b and x >= 0.
 *
 * A is row-major with constraint_count rows and variable_count columns.
 * A and b may be NULL only when constraint_count is zero. The solution and
 * objective outputs are written only for an optimal result. iterations is
 * optional. Equalities and finite upper bounds can be encoded as pairs of
 * inequalities.
 */
numelyx_core_lp_status numelyx_core_lp_default_options(
    numelyx_core_lp_options *options
);

numelyx_core_lp_status numelyx_core_lp_solve_maximize(
    const double *a,
    const double *b,
    const double *c,
    size_t constraint_count,
    size_t variable_count,
    const numelyx_core_lp_options *options,
    double *solution,
    double *objective_value,
    size_t *iterations
);

#endif
