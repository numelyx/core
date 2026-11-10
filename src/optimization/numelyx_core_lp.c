/**
 * @file numelyx_core_lp.c
 * @brief Two-phase tableau simplex with Bland anti-cycling selection.
 */
#include "numelyx_core_lp.h"
#include "numelyx_core_size.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    SIMPLEX_OPTIMAL = 0,
    SIMPLEX_UNBOUNDED,
    SIMPLEX_LIMIT,
    SIMPLEX_NUMERICAL_FAILURE
} simplex_status;

typedef struct {
    size_t constraints;
    size_t variables;
    size_t rows;
    size_t columns;
    size_t *basis;
    size_t *nonbasis;
    double *tableau;
    const numelyx_core_lp_options *options;
    size_t iterations;
} simplex_workspace;

static double *tableau_cell(simplex_workspace *workspace, size_t row, size_t col)
{
    return &workspace->tableau[row * workspace->columns + col];
}

static bool finite_tableau(const simplex_workspace *workspace)
{
    size_t index;
    size_t count = workspace->rows * workspace->columns;
    for (index = 0; index < count; index++) {
        if (!isfinite(workspace->tableau[index])) return false;
    }
    return true;
}

static bool simplex_pivot(simplex_workspace *workspace, size_t row, size_t col)
{
    size_t i;
    size_t j;
    double pivot = *tableau_cell(workspace, row, col);
    double inverse;

    if (!isfinite(pivot) || pivot == 0.0) return false;
    inverse = 1.0 / pivot;
    if (!isfinite(inverse)) return false;

    for (i = 0; i < workspace->rows; i++) {
        if (i == row) continue;
        for (j = 0; j < workspace->columns; j++) {
            double *cell;
            double value;
            if (j == col) continue;
            cell = tableau_cell(workspace, i, j);
            value = *cell -
                (*tableau_cell(workspace, row, j)) *
                (*tableau_cell(workspace, i, col)) * inverse;
            if (!isfinite(value)) return false;
            *cell = value;
        }
    }

    for (j = 0; j < workspace->columns; j++) {
        if (j == col) continue;
        *tableau_cell(workspace, row, j) *= inverse;
        if (!isfinite(*tableau_cell(workspace, row, j))) return false;
    }
    for (i = 0; i < workspace->rows; i++) {
        if (i == row) continue;
        *tableau_cell(workspace, i, col) *= -inverse;
        if (!isfinite(*tableau_cell(workspace, i, col))) return false;
    }
    *tableau_cell(workspace, row, col) = inverse;

    {
        size_t temporary = workspace->basis[row];
        workspace->basis[row] = workspace->nonbasis[col];
        workspace->nonbasis[col] = temporary;
    }
    return finite_tableau(workspace);
}

static simplex_status simplex_optimize(
    simplex_workspace *workspace, bool phase_one
)
{
    size_t objective_row = phase_one
        ? workspace->constraints + 1
        : workspace->constraints;

    for (;;) {
        size_t entering = (size_t) -1;
        size_t leaving = (size_t) -1;
        size_t col;
        size_t row;
        double best_ratio = 0.0;

        /* Bland's rule: choose the smallest-index eligible nonbasic variable. */
        for (col = 0; col <= workspace->variables; col++) {
            if (!phase_one && workspace->nonbasis[col] == (size_t) -1)
                continue;
            if (*tableau_cell(workspace, objective_row, col) <
                -workspace->options->optimality_tolerance) {
                if (entering == (size_t) -1 ||
                    workspace->nonbasis[col] < workspace->nonbasis[entering]) {
                    entering = col;
                }
            }
        }

        if (entering == (size_t) -1) return SIMPLEX_OPTIMAL;
        if (workspace->iterations >= workspace->options->max_iterations)
            return SIMPLEX_LIMIT;

        for (row = 0; row < workspace->constraints; row++) {
            double coefficient = *tableau_cell(workspace, row, entering);
            double rhs;
            double ratio;
            double tolerance;

            if (coefficient <= workspace->options->feasibility_tolerance)
                continue;
            rhs = *tableau_cell(workspace, row, workspace->variables + 1);
            if (rhs < -workspace->options->feasibility_tolerance)
                return SIMPLEX_NUMERICAL_FAILURE;
            if (rhs < 0.0) rhs = 0.0;
            ratio = rhs / coefficient;
            if (!isfinite(ratio)) return SIMPLEX_NUMERICAL_FAILURE;

            if (leaving == (size_t) -1) {
                leaving = row;
                best_ratio = ratio;
                continue;
            }
            tolerance = workspace->options->feasibility_tolerance *
                fmax(1.0, fmax(fabs(best_ratio), fabs(ratio)));
            if (ratio < best_ratio - tolerance ||
                (fabs(ratio - best_ratio) <= tolerance &&
                 workspace->basis[row] < workspace->basis[leaving])) {
                leaving = row;
                best_ratio = ratio;
            }
        }

        if (leaving == (size_t) -1) return SIMPLEX_UNBOUNDED;
        if (!simplex_pivot(workspace, leaving, entering))
            return SIMPLEX_NUMERICAL_FAILURE;
        workspace->iterations++;
    }
}

numelyx_core_lp_status numelyx_core_lp_default_options(
    numelyx_core_lp_options *options
)
{
    if (options == NULL) return NUMELYX_LP_INVALID_ARGUMENT;
    options->max_iterations = 10000;
    options->feasibility_tolerance = 1e-9;
    options->optimality_tolerance = 1e-9;
    return NUMELYX_LP_SUCCESS;
}

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
)
{
    simplex_workspace workspace;
    size_t tableau_rows;
    size_t tableau_columns;
    size_t tableau_elements;
    size_t tableau_bytes;
    size_t basis_bytes = 0;
    size_t nonbasis_bytes;
    size_t candidate_bytes;
    size_t row;
    size_t col;
    double *candidate = NULL;
    double *tableau = NULL;
    size_t *basis = NULL;
    size_t *nonbasis = NULL;
    numelyx_core_lp_status result_status = NUMELYX_LP_NUMERICAL_FAILURE;
    simplex_status phase_status;

    if (solution == NULL || objective_value == NULL || c == NULL ||
        variable_count == 0 || options == NULL ||
        options->max_iterations == 0) {
        return NUMELYX_LP_INVALID_ARGUMENT;
    }
    if (!isfinite(options->feasibility_tolerance) ||
        !isfinite(options->optimality_tolerance)) {
        return NUMELYX_LP_NON_FINITE_INPUT;
    }
    if (options->feasibility_tolerance <= 0.0 ||
        options->optimality_tolerance <= 0.0) {
        return NUMELYX_LP_INVALID_ARGUMENT;
    }
    if (constraint_count > 0 && (a == NULL || b == NULL))
        return NUMELYX_LP_INVALID_ARGUMENT;

    for (col = 0; col < variable_count; col++) {
        if (!isfinite(c[col])) return NUMELYX_LP_NON_FINITE_INPUT;
    }
    for (row = 0; row < constraint_count; row++) {
        if (!isfinite(b[row])) return NUMELYX_LP_NON_FINITE_INPUT;
        for (col = 0; col < variable_count; col++) {
            if (!isfinite(a[row * variable_count + col]))
                return NUMELYX_LP_NON_FINITE_INPUT;
        }
    }

    if (constraint_count > SIZE_MAX - 2 ||
        variable_count > SIZE_MAX - 2) {
        return NUMELYX_LP_SIZE_OVERFLOW;
    }
    tableau_rows = constraint_count + 2;
    tableau_columns = variable_count + 2;
    if (!numelyx_core_size_multiply(
            tableau_rows, tableau_columns, &tableau_elements) ||
        !numelyx_core_size_multiply(
            tableau_elements, sizeof(*tableau), &tableau_bytes) ||
        !numelyx_core_size_multiply(
            variable_count + 1, sizeof(*nonbasis), &nonbasis_bytes) ||
        !numelyx_core_size_multiply(
            variable_count, sizeof(*candidate), &candidate_bytes) ||
        (constraint_count > 0 &&
         !numelyx_core_size_multiply(
             constraint_count, sizeof(*basis), &basis_bytes))) {
        return NUMELYX_LP_SIZE_OVERFLOW;
    }

    tableau = (double *) calloc(1, tableau_bytes);
    nonbasis = (size_t *) malloc(nonbasis_bytes);
    candidate = (double *) malloc(candidate_bytes);
    if (constraint_count > 0)
        basis = (size_t *) malloc(basis_bytes);
    if (tableau == NULL || nonbasis == NULL || candidate == NULL ||
        (constraint_count > 0 && basis == NULL)) {
        result_status = NUMELYX_LP_OUT_OF_MEMORY;
        goto cleanup;
    }

    workspace.constraints = constraint_count;
    workspace.variables = variable_count;
    workspace.rows = tableau_rows;
    workspace.columns = tableau_columns;
    workspace.basis = basis;
    workspace.nonbasis = nonbasis;
    workspace.tableau = tableau;
    workspace.options = options;
    workspace.iterations = 0;

    for (row = 0; row < constraint_count; row++) {
        for (col = 0; col < variable_count; col++) {
            *tableau_cell(&workspace, row, col) =
                a[row * variable_count + col];
        }
        workspace.basis[row] = variable_count + row;
        *tableau_cell(&workspace, row, variable_count) = -1.0;
        *tableau_cell(&workspace, row, variable_count + 1) = b[row];
    }
    for (col = 0; col < variable_count; col++) {
        workspace.nonbasis[col] = col;
        *tableau_cell(&workspace, constraint_count, col) = -c[col];
    }
    workspace.nonbasis[variable_count] = (size_t) -1;
    *tableau_cell(&workspace, constraint_count + 1, variable_count) = 1.0;

    if (constraint_count > 0) {
        size_t most_negative_row = 0;
        for (row = 1; row < constraint_count; row++) {
            if (*tableau_cell(&workspace, row, variable_count + 1) <
                *tableau_cell(&workspace, most_negative_row, variable_count + 1)) {
                most_negative_row = row;
            }
        }

        if (*tableau_cell(&workspace, most_negative_row, variable_count + 1) <
            -options->feasibility_tolerance) {
            if (!simplex_pivot(&workspace, most_negative_row, variable_count)) {
                result_status = NUMELYX_LP_NUMERICAL_FAILURE;
                goto cleanup;
            }
            phase_status = simplex_optimize(&workspace, true);
            if (phase_status == SIMPLEX_LIMIT) {
                result_status = NUMELYX_LP_ITERATION_LIMIT;
                goto cleanup;
            }
            if (phase_status != SIMPLEX_OPTIMAL) {
                result_status = NUMELYX_LP_NUMERICAL_FAILURE;
                goto cleanup;
            }
            if (*tableau_cell(
                    &workspace, constraint_count + 1, variable_count + 1
                ) < -options->feasibility_tolerance) {
                result_status = NUMELYX_LP_INFEASIBLE;
                goto cleanup;
            }
            if (fabs(*tableau_cell(
                    &workspace, constraint_count + 1, variable_count + 1
                )) > options->feasibility_tolerance) {
                result_status = NUMELYX_LP_NUMERICAL_FAILURE;
                goto cleanup;
            }

            for (row = 0; row < constraint_count; row++) {
                if (workspace.basis[row] == (size_t) -1) {
                    size_t entering = (size_t) -1;
                    double largest = 0.0;
                    for (col = 0; col < variable_count; col++) {
                        double magnitude = fabs(*tableau_cell(&workspace, row, col));
                        if (magnitude > largest) {
                            largest = magnitude;
                            entering = col;
                        }
                    }
                    if (entering != (size_t) -1 &&
                        largest > options->feasibility_tolerance) {
                        if (!simplex_pivot(&workspace, row, entering)) {
                            result_status = NUMELYX_LP_NUMERICAL_FAILURE;
                            goto cleanup;
                        }
                    }
                }
            }
        }
    }

    phase_status = simplex_optimize(&workspace, false);
    if (phase_status == SIMPLEX_UNBOUNDED) {
        result_status = NUMELYX_LP_UNBOUNDED;
        goto cleanup;
    }
    if (phase_status == SIMPLEX_LIMIT) {
        result_status = NUMELYX_LP_ITERATION_LIMIT;
        goto cleanup;
    }
    if (phase_status != SIMPLEX_OPTIMAL) {
        result_status = NUMELYX_LP_NUMERICAL_FAILURE;
        goto cleanup;
    }

    memset(candidate, 0, candidate_bytes);
    for (row = 0; row < constraint_count; row++) {
        if (workspace.basis[row] < variable_count) {
            double value = *tableau_cell(&workspace, row, variable_count + 1);
            if (value < -options->feasibility_tolerance || !isfinite(value)) {
                result_status = NUMELYX_LP_NUMERICAL_FAILURE;
                goto cleanup;
            }
            candidate[workspace.basis[row]] = value < 0.0 ? 0.0 : value;
        }
    }
    {
        double objective = *tableau_cell(
            &workspace, constraint_count, variable_count + 1);
        if (!isfinite(objective)) {
            result_status = NUMELYX_LP_NUMERICAL_FAILURE;
            goto cleanup;
        }
        memcpy(solution, candidate, candidate_bytes);
        *objective_value = objective;
    }
    if (iterations != NULL) *iterations = workspace.iterations;
    result_status = NUMELYX_LP_OPTIMAL;

cleanup:
    free(tableau);
    free(nonbasis);
    free(candidate);
    free(basis);
    return result_status;
}
