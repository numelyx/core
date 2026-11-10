/**
 * @file numelyx_core_matrix_condition.c
 * @brief Matrix condition-number estimation using reusable LU solves.
 *
 * @details The inverse 1-norm is estimated with a Hager-style iteration.
 * The implementation reuses one cached LU factorization and never forms the
 * inverse matrix.
 */

#include "numelyx_core_matrix.h"
#include "numelyx_core_matrix_internal.h"
#include "numelyx_core_matrix_lu.h"
#include "numelyx_core_size.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static numelyx_core_matrix_status condition_scaled_one_norm(
    const numelyx_core_matrix *matrix,
    size_t size,
    double *matrix_scale,
    double *scaled_norm
)
{
    size_t row;
    size_t column;
    double scale = 0.0;
    double norm = 0.0;

    for (row = 0; row < size; row++) {
        for (column = 0; column < size; column++) {
            double value;
            numelyx_core_matrix_status status = numelyx_core_matrix_get(
                matrix, row, column, &value
            );

            if (status != NUMERUS_MATRIX_SUCCESS) {
                return status;
            }
            if (!isfinite(value)) {
                return NUMERUS_MATRIX_NON_FINITE;
            }
            if (fabs(value) > scale) {
                scale = fabs(value);
            }
        }
    }

    if (scale == 0.0) {
        return NUMERUS_MATRIX_SINGULAR;
    }

    for (column = 0; column < size; column++) {
        double column_sum = 0.0;

        for (row = 0; row < size; row++) {
            double value;
            numelyx_core_matrix_status status = numelyx_core_matrix_get(
                matrix, row, column, &value
            );

            if (status != NUMERUS_MATRIX_SUCCESS) {
                return status;
            }

            column_sum += fabs(value / scale);
        }

        if (column_sum > norm) {
            norm = column_sum;
        }
    }

    *matrix_scale = scale;
    *scaled_norm = norm;
    return NUMERUS_MATRIX_SUCCESS;
}

static numelyx_core_matrix_status condition_solve(
    const numelyx_core_matrix_lu_factorization *factorization,
    size_t size,
    double matrix_scale,
    const double *right_hand_side,
    double *solution,
    double *work
)
{
    size_t row;

    /* Forward substitution: L y = P b. */
    for (row = 0; row < size; row++) {
        size_t source_row;
        size_t column;
        double value;

        {
            numelyx_core_matrix_status status =
                numelyx_core_matrix_lu_get_permutation(
                    factorization, row, &source_row
                );

            if (status != NUMERUS_MATRIX_SUCCESS) {
                return status;
            }
        }

        value = matrix_scale * right_hand_side[source_row];
        if (!isfinite(value)) {
            return NUMERUS_MATRIX_NON_FINITE;
        }

        for (column = 0; column < row; column++) {
            double lower_value;
            numelyx_core_matrix_status status =
                numelyx_core_matrix_lu_get_lower(
                    factorization, row, column, &lower_value
                );

            if (status != NUMERUS_MATRIX_SUCCESS) {
                return status;
            }

            value -= lower_value * work[column];
            if (!isfinite(value)) {
                return NUMERUS_MATRIX_NON_FINITE;
            }
        }

        work[row] = value;
    }

    /* Back substitution: U x = y. */
    row = size;
    while (row > 0) {
        size_t column;
        double value;
        double diagonal;
        numelyx_core_matrix_status status;

        row--;
        value = work[row];

        for (column = row + 1; column < size; column++) {
            double upper_value;

            status = numelyx_core_matrix_lu_get_upper(
                factorization, row, column, &upper_value
            );
            if (status != NUMERUS_MATRIX_SUCCESS) {
                return status;
            }

            value -= upper_value * solution[column];
            if (!isfinite(value)) {
                return NUMERUS_MATRIX_NON_FINITE;
            }
        }

        status = numelyx_core_matrix_lu_get_upper(
            factorization, row, row, &diagonal
        );
        if (status != NUMERUS_MATRIX_SUCCESS) {
            return status;
        }
        if (diagonal == 0.0) {
            return NUMERUS_MATRIX_SINGULAR;
        }
        if (!isfinite(diagonal)) {
            return NUMERUS_MATRIX_NON_FINITE;
        }

        value /= diagonal;
        if (!isfinite(value)) {
            return NUMERUS_MATRIX_NON_FINITE;
        }

        solution[row] = value;
    }

    return NUMERUS_MATRIX_SUCCESS;
}

static numelyx_core_matrix_status condition_solve_transpose(
    const numelyx_core_matrix_lu_factorization *factorization,
    size_t size,
    double matrix_scale,
    const double *right_hand_side,
    double *solution,
    double *work
)
{
    size_t row;

    /* Forward substitution: U^T y = b. */
    for (row = 0; row < size; row++) {
        size_t column;
        double value = matrix_scale * right_hand_side[row];
        double diagonal;
        numelyx_core_matrix_status status;

        if (!isfinite(value)) {
            return NUMERUS_MATRIX_NON_FINITE;
        }

        for (column = 0; column < row; column++) {
            double upper_value;

            status = numelyx_core_matrix_lu_get_upper(
                factorization, column, row, &upper_value
            );
            if (status != NUMERUS_MATRIX_SUCCESS) {
                return status;
            }

            value -= upper_value * work[column];
            if (!isfinite(value)) {
                return NUMERUS_MATRIX_NON_FINITE;
            }
        }

        status = numelyx_core_matrix_lu_get_upper(
            factorization, row, row, &diagonal
        );
        if (status != NUMERUS_MATRIX_SUCCESS) {
            return status;
        }
        if (diagonal == 0.0) {
            return NUMERUS_MATRIX_SINGULAR;
        }
        if (!isfinite(diagonal)) {
            return NUMERUS_MATRIX_NON_FINITE;
        }

        value /= diagonal;
        if (!isfinite(value)) {
            return NUMERUS_MATRIX_NON_FINITE;
        }

        work[row] = value;
    }

    /* Back substitution: L^T z = y (L has a unit diagonal). */
    row = size;
    while (row > 0) {
        size_t column;
        double value;

        row--;
        value = work[row];

        for (column = row + 1; column < size; column++) {
            double lower_value;
            numelyx_core_matrix_status status =
                numelyx_core_matrix_lu_get_lower(
                    factorization, column, row, &lower_value
                );

            if (status != NUMERUS_MATRIX_SUCCESS) {
                return status;
            }

            value -= lower_value * solution[column];
            if (!isfinite(value)) {
                return NUMERUS_MATRIX_NON_FINITE;
            }
        }

        solution[row] = value;
    }

    /* P x = z, so x[permutation[row]] = z[row]. */
    for (row = 0; row < size; row++) {
        size_t source_row;
        numelyx_core_matrix_status status =
            numelyx_core_matrix_lu_get_permutation(
                factorization, row, &source_row
            );

        if (status != NUMERUS_MATRIX_SUCCESS) {
            return status;
        }
        work[row] = solution[row];
    }

    for (row = 0; row < size; row++) {
        size_t source_row;
        numelyx_core_matrix_status status =
            numelyx_core_matrix_lu_get_permutation(
                factorization, row, &source_row
            );

        if (status != NUMERUS_MATRIX_SUCCESS) {
            return status;
        }
        solution[source_row] = work[row];
    }

    return NUMERUS_MATRIX_SUCCESS;
}

numelyx_core_matrix_status numelyx_core_matrix_condition_estimate_one(
    const numelyx_core_matrix *matrix,
    double *condition_estimate
)
{
    numelyx_core_matrix_lu_factorization *factorization = NULL;
    size_t size;
    size_t rank;
    size_t workspace_count;
    size_t workspace_bytes;
    size_t index;
    size_t iteration;
    size_t last_index = SIZE_MAX;
    double matrix_scale;
    double scaled_matrix_norm;
    double inverse_norm = 0.0;
    double *workspace = NULL;
    double *x;
    double *y;
    double *signs;
    double *z;
    double *work;
    numelyx_core_matrix_status status;

    if (matrix == NULL || condition_estimate == NULL) {
        return NUMERUS_MATRIX_INVALID_ARGUMENT;
    }
    if (numelyx_core_matrix_rows(matrix) != numelyx_core_matrix_columns(matrix)) {
        return NUMERUS_MATRIX_NOT_SQUARE;
    }

    size = numelyx_core_matrix_rows(matrix);
    if (size == 0) {
        return NUMERUS_MATRIX_INVALID_ARGUMENT;
    }

    status = numelyx_core_matrix_rank(matrix, &rank);
    if (status != NUMERUS_MATRIX_SUCCESS) {
        return status;
    }
    if (rank != size) {
        *condition_estimate = INFINITY;
        return NUMERUS_MATRIX_SUCCESS;
    }

    if (!numelyx_core_size_multiply(size, 5, &workspace_count) ||
        !numelyx_core_size_multiply(
            workspace_count, sizeof(*workspace), &workspace_bytes
        )) {
        return NUMERUS_MATRIX_OVERFLOW;
    }

    workspace = malloc(workspace_bytes);
    if (workspace == NULL) {
        return NUMERUS_MATRIX_OUT_OF_MEMORY;
    }

    x = workspace;
    y = x + size;
    signs = y + size;
    z = signs + size;
    work = z + size;

    status = numelyx_core_matrix_get_or_factorize_lu(
        matrix, &factorization
    );
    if (status != NUMERUS_MATRIX_SUCCESS) {
        goto cleanup;
    }
    if (numelyx_core_matrix_lu_rank(factorization) != size) {
        *condition_estimate = INFINITY;
        status = NUMERUS_MATRIX_SUCCESS;
        goto cleanup;
    }

    status = condition_scaled_one_norm(
        matrix, size, &matrix_scale, &scaled_matrix_norm
    );
    if (status != NUMERUS_MATRIX_SUCCESS) {
        goto cleanup;
    }

    {
        double initial_component = 1.0 / (double) size;

        if (matrix_scale * initial_component == 0.0) {
            for (index = 0; index < size; index++) {
                x[index] = 0.0;
            }
            x[0] = 1.0;
        } else {
            for (index = 0; index < size; index++) {
                x[index] = initial_component;
            }
        }
    }

    /*
     * Hager's estimator alternates solves with A and A^T. The iteration is
     * capped at five passes to bound overhead; as a norm estimator, its result
     * may be lower than the exact inverse 1-norm. Each pass costs O(n^2) after
     * the cached LU factorization.
     */
    for (iteration = 0; iteration < 5; iteration++) {
        double estimate = 0.0;
        double maximum = -1.0;
        double dual_product = 0.0;
        size_t maximum_index = 0;

        status = condition_solve(
            factorization, size, matrix_scale, x, y, work
        );
        if (status != NUMERUS_MATRIX_SUCCESS) {
            goto cleanup;
        }

        for (index = 0; index < size; index++) {
            estimate += fabs(y[index]);
            signs[index] = y[index] >= 0.0 ? 1.0 : -1.0;
        }
        if (!isfinite(estimate)) {
            inverse_norm = INFINITY;
            break;
        }
        if (estimate > inverse_norm) {
            inverse_norm = estimate;
        }

        status = condition_solve_transpose(
            factorization, size, matrix_scale, signs, z, work
        );
        if (status != NUMERUS_MATRIX_SUCCESS) {
            goto cleanup;
        }

        for (index = 0; index < size; index++) {
            double magnitude = fabs(z[index]);

            dual_product += z[index] * x[index];
            if (magnitude > maximum) {
                maximum = magnitude;
                maximum_index = index;
            }
        }
        if (!isfinite(dual_product)) {
            inverse_norm = INFINITY;
            break;
        }
        if (maximum <= dual_product || maximum_index == last_index) {
            break;
        }

        for (index = 0; index < size; index++) {
            x[index] = 0.0;
        }
        x[maximum_index] = 1.0;
        last_index = maximum_index;
    }

    *condition_estimate = scaled_matrix_norm * inverse_norm;
    status = NUMERUS_MATRIX_SUCCESS;

cleanup:
    if (factorization != NULL) {
        numelyx_core_matrix_lu_destroy(factorization);
    }
    free(workspace);
    return status;
}
