#include <numelyx/matrix.h>

#include <stdint.h>
#include <stdlib.h>

#include "numelyx_core_matrix_core.h"
#include "numelyx_core_matrix_types.h"
#include "numelyx_core_matrix_views.h"
#include "numelyx_core_matrix_operations.h"
#include "numelyx_core_matrix_linalg.h"

static numelyx_status map_matrix_status(numelyx_core_matrix_status status)
{
    switch (status) {
        case NUMERUS_MATRIX_SUCCESS:
            return NUMELYX_SUCCESS;
        case NUMERUS_MATRIX_INVALID_ARGUMENT:
            return NUMELYX_INVALID_ARGUMENT;
        case NUMERUS_MATRIX_OUT_OF_MEMORY:
            return NUMELYX_OUT_OF_MEMORY;
        case NUMERUS_MATRIX_OVERFLOW:
            return NUMELYX_SIZE_OVERFLOW;
        case NUMERUS_MATRIX_OUT_OF_BOUNDS:
            return NUMELYX_OUT_OF_BOUNDS;
        case NUMERUS_MATRIX_NOT_SQUARE:
        case NUMERUS_MATRIX_DIMENSION_MISMATCH:
            return NUMELYX_DIMENSION_MISMATCH;
        case NUMERUS_MATRIX_NON_FINITE:
            return NUMELYX_NON_FINITE;
        case NUMERUS_MATRIX_SINGULAR:
            return NUMELYX_SINGULAR;
        case NUMERUS_MATRIX_RANK_DEFICIENT:
            return NUMELYX_RANK_DEFICIENT;
        case NUMERUS_MATRIX_NO_CONVERGENCE:
            return NUMELYX_NO_CONVERGENCE;
        case NUMERUS_MATRIX_DIVISION_BY_ZERO:
            return NUMELYX_DIVISION_BY_ZERO;
        case NUMERUS_MATRIX_NOT_SYMMETRIC:
            return NUMELYX_NOT_SYMMETRIC;
        case NUMERUS_MATRIX_NOT_POSITIVE_DEFINITE:
            return NUMELYX_NOT_POSITIVE_DEFINITE;
        case NUMERUS_MATRIX_NOT_POSITIVE_SEMIDEFINITE:
            return NUMELYX_NOT_POSITIVE_SEMIDEFINITE;
        case NUMERUS_MATRIX_PIVOT_TOO_SMALL:
            return NUMELYX_PIVOT_TOO_SMALL;
        default:
            return NUMELYX_NUMERICAL_ERROR;
    }
}

static numelyx_status create_matrix(
    numelyx_core_storage_kind kind,
    size_t rows,
    size_t columns,
    const numelyx_core_matrix_data *data,
    numelyx_matrix **matrix
)
{
    int status;

    if (matrix == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }

    *matrix = NULL;
    status = numelyx_core_matrix_create(kind, rows, columns, data, matrix);
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_create_dense(
    size_t rows,
    size_t columns,
    const double *values,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};

    if (values == NULL) {
        if (matrix != NULL) {
            *matrix = NULL;
        }
        return NUMELYX_INVALID_ARGUMENT;
    }

    data.values = values;
    return create_matrix(NUMERUS_STORAGE_DENSE, rows, columns, &data, matrix);
}

numelyx_status numelyx_matrix_create_identity(
    size_t size,
    numelyx_matrix **matrix
)
{
    return create_matrix(NUMERUS_STORAGE_IDENTITY, size, size, NULL, matrix);
}

numelyx_status numelyx_matrix_create_zero(
    size_t rows,
    size_t columns,
    numelyx_matrix **matrix
)
{
    return create_matrix(NUMERUS_STORAGE_ZERO, rows, columns, NULL, matrix);
}

numelyx_status numelyx_matrix_create_constant(
    size_t rows,
    size_t columns,
    double value,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};

    data.value = value;
    return create_matrix(NUMERUS_STORAGE_CONSTANT, rows, columns, &data, matrix);
}

numelyx_status numelyx_matrix_create_diagonal(
    size_t size,
    const double *diagonal,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};

    if (diagonal == NULL) {
        if (matrix != NULL) {
            *matrix = NULL;
        }
        return NUMELYX_INVALID_ARGUMENT;
    }

    data.values = diagonal;
    return create_matrix(NUMERUS_STORAGE_DIAGONAL, size, size, &data, matrix);
}

numelyx_status numelyx_matrix_create_scaled_identity(
    size_t size,
    double diagonal_value,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};

    data.value = diagonal_value;
    return create_matrix(
        NUMERUS_STORAGE_SCALED_IDENTITY,
        size,
        size,
        &data,
        matrix
    );
}

numelyx_status numelyx_matrix_create_upper_triangular(
    size_t size,
    const double *packed_values,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};

    if (packed_values == NULL) {
        if (matrix != NULL) {
            *matrix = NULL;
        }
        return NUMELYX_INVALID_ARGUMENT;
    }
    data.values = packed_values;
    return create_matrix(
        NUMERUS_STORAGE_UPPER_TRIANGULAR,
        size,
        size,
        &data,
        matrix
    );
}

numelyx_status numelyx_matrix_create_lower_triangular(
    size_t size,
    const double *packed_values,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};

    if (packed_values == NULL) {
        if (matrix != NULL) {
            *matrix = NULL;
        }
        return NUMELYX_INVALID_ARGUMENT;
    }
    data.values = packed_values;
    return create_matrix(
        NUMERUS_STORAGE_LOWER_TRIANGULAR,
        size,
        size,
        &data,
        matrix
    );
}

numelyx_status numelyx_matrix_create_symmetric(
    size_t size,
    const double *packed_lower_values,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};

    if (packed_lower_values == NULL) {
        if (matrix != NULL) {
            *matrix = NULL;
        }
        return NUMELYX_INVALID_ARGUMENT;
    }
    data.values = packed_lower_values;
    return create_matrix(
        NUMERUS_STORAGE_SYMMETRIC,
        size,
        size,
        &data,
        matrix
    );
}

numelyx_status numelyx_matrix_create_banded(
    size_t rows,
    size_t columns,
    size_t lower_bandwidth,
    size_t upper_bandwidth,
    const double *packed_band_values,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};

    if (matrix == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *matrix = NULL;

    if (rows == 0 || columns == 0
        || lower_bandwidth >= rows
        || upper_bandwidth >= columns
        || packed_band_values == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }

    data.values = packed_band_values;
    data.lower_bandwidth = lower_bandwidth;
    data.upper_bandwidth = upper_bandwidth;
    return create_matrix(
        NUMERUS_STORAGE_BANDED,
        rows,
        columns,
        &data,
        matrix
    );
}

numelyx_status numelyx_matrix_create_sparse(
    size_t rows,
    size_t columns,
    double default_value,
    const numelyx_sparse_entry *entries,
    size_t count,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};
    numelyx_core_storage_sparse_entry *core_entries = NULL;
    numelyx_status result;

    if (matrix == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *matrix = NULL;

    if (count != 0 && entries == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    if (count > SIZE_MAX / sizeof(*core_entries)) {
        return NUMELYX_SIZE_OVERFLOW;
    }

    if (count != 0) {
        core_entries = malloc(count * sizeof(*core_entries));
        if (core_entries == NULL) {
            return NUMELYX_OUT_OF_MEMORY;
        }

        for (size_t i = 0; i < count; i++) {
            core_entries[i].index = entries[i].index;
            core_entries[i].value = entries[i].value;
        }
    }

    data.entries = core_entries;
    data.count = count;
    data.default_value = default_value;
    result = create_matrix(
        NUMERUS_STORAGE_SPARSE,
        rows,
        columns,
        &data,
        matrix
    );

    free(core_entries);
    return result;
}

numelyx_status numelyx_matrix_transpose(
    numelyx_matrix *parent,
    numelyx_matrix **view
)
{
    int status;

    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;

    status = numelyx_core_matrix_create_transpose(parent, view);
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_slice(
    numelyx_matrix *parent,
    size_t row_start,
    size_t row_count,
    size_t column_start,
    size_t column_count,
    numelyx_matrix **view
)
{
    int status;

    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;

    status = numelyx_core_matrix_create_slice(
        parent,
        row_start,
        row_count,
        column_start,
        column_count,
        view
    );
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_flip_rows(
    numelyx_matrix *parent,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_flip_rows(parent, view));
}

numelyx_status numelyx_matrix_flip_columns(
    numelyx_matrix *parent,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_flip_columns(parent, view));
}

numelyx_status numelyx_matrix_select_rows(
    numelyx_matrix *parent,
    const size_t *indices,
    size_t count,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_select_rows(parent, indices, count, view));
}

numelyx_status numelyx_matrix_select_columns(
    numelyx_matrix *parent,
    const size_t *indices,
    size_t count,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_select_columns(parent, indices, count, view));
}

numelyx_status numelyx_matrix_reshape(
    numelyx_matrix *parent,
    size_t rows,
    size_t columns,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_reshape(parent, rows, columns, view));
}

numelyx_status numelyx_matrix_flatten(
    numelyx_matrix *parent,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_flatten(parent, view));
}

numelyx_status numelyx_matrix_remove_row(
    numelyx_matrix *parent,
    size_t row,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_remove_row(parent, row, view));
}

numelyx_status numelyx_matrix_remove_column(
    numelyx_matrix *parent,
    size_t column,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_remove_column(parent, column, view));
}

numelyx_status numelyx_matrix_swap_rows(
    numelyx_matrix *parent,
    size_t row1,
    size_t row2,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_swap_rows(parent, row1, row2, view));
}

numelyx_status numelyx_matrix_swap_columns(
    numelyx_matrix *parent,
    size_t column1,
    size_t column2,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_swap_columns(parent, column1, column2, view));
}

numelyx_status numelyx_matrix_rotate_90_clockwise(
    numelyx_matrix *parent,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_rotate_90_clockwise(parent, view));
}

numelyx_status numelyx_matrix_rotate_180(
    numelyx_matrix *parent,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_rotate_180(parent, view));
}

numelyx_status numelyx_matrix_rotate_90_counterclockwise(
    numelyx_matrix *parent,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_rotate_90_counterclockwise(parent, view));
}

numelyx_status numelyx_matrix_join_horizontal(
    numelyx_matrix *left,
    numelyx_matrix *right,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_join_horizontal(left, right, view));
}

numelyx_status numelyx_matrix_join_vertical(
    numelyx_matrix *top,
    numelyx_matrix *bottom,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_join_vertical(top, bottom, view));
}

numelyx_status numelyx_matrix_pad(
    numelyx_matrix *parent,
    size_t top,
    size_t bottom,
    size_t left,
    size_t right,
    double value,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_pad(parent, top, bottom, left, right, value, view));
}

numelyx_status numelyx_matrix_zero_extend(
    numelyx_matrix *parent,
    size_t top,
    size_t bottom,
    size_t left,
    size_t right,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_zero_extend(parent, top, bottom, left, right, view));
}

numelyx_status numelyx_matrix_repeat(
    numelyx_matrix *parent,
    size_t row_repetitions,
    size_t column_repetitions,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_repeat(parent, row_repetitions, column_repetitions, view));
}

numelyx_status numelyx_matrix_block_diagonal(
    numelyx_matrix *parent,
    size_t repetitions,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_block_diagonal(parent, repetitions, view));
}

numelyx_status numelyx_matrix_block_grid(
    numelyx_matrix *const *blocks,
    size_t block_row_count,
    size_t block_column_count,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_block_grid(blocks, block_row_count, block_column_count, view));
}

numelyx_status numelyx_matrix_diagonal_extract(
    numelyx_matrix *parent,
    ptrdiff_t offset,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_diagonal_extract(parent, offset, view));
}

numelyx_status numelyx_matrix_diagonal_from_vector(
    numelyx_matrix *vector,
    ptrdiff_t offset,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_diagonal_from_vector(vector, offset, view));
}

numelyx_status numelyx_matrix_permute_rows(
    numelyx_matrix *parent,
    const size_t *permutation,
    size_t count,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_permute_rows(parent, permutation, count, view));
}

numelyx_status numelyx_matrix_permute_columns(
    numelyx_matrix *parent,
    const size_t *permutation,
    size_t count,
    numelyx_matrix **view
)
{
    if (view == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *view = NULL;
    return map_matrix_status((numelyx_core_matrix_status)
        numelyx_core_matrix_create_permute_columns(parent, permutation, count, view));
}

numelyx_status numelyx_matrix_add(
    numelyx_matrix *left,
    numelyx_matrix *right,
    numelyx_matrix **result
)
{
    int status;

    if (result == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *result = NULL;
    status = numelyx_core_matrix_create_add(left, right, result);
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_subtract(
    numelyx_matrix *left,
    numelyx_matrix *right,
    numelyx_matrix **result
)
{
    int status;

    if (result == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *result = NULL;
    status = numelyx_core_matrix_create_subtract(left, right, result);
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_hadamard_product(
    numelyx_matrix *left,
    numelyx_matrix *right,
    numelyx_matrix **result
)
{
    int status;

    if (result == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *result = NULL;
    status = numelyx_core_matrix_create_hadamard_product(left, right, result);
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_multiply(
    const numelyx_matrix *left,
    const numelyx_matrix *right,
    numelyx_matrix **result
)
{
    int status;

    if (result == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *result = NULL;
    status = numelyx_core_matrix_multiply(left, right, result);
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_scale(
    numelyx_matrix *source,
    double scalar,
    numelyx_matrix **result
)
{
    int status;

    if (result == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *result = NULL;
    status = numelyx_core_matrix_create_scale(source, scalar, result);
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_materialize(
    const numelyx_matrix *source,
    numelyx_matrix **result
)
{
    int status;

    if (result == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *result = NULL;
    status = numelyx_core_matrix_materialize(source, result);
    return map_matrix_status((numelyx_core_matrix_status) status);
}

numelyx_status numelyx_matrix_determinant(
    numelyx_matrix *matrix,
    double *determinant
)
{
    return map_matrix_status(numelyx_core_matrix_determinant(
        matrix,
        determinant
    ));
}

numelyx_status numelyx_matrix_inverse(
    const numelyx_matrix *matrix,
    numelyx_matrix **inverse
)
{
    if (inverse == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *inverse = NULL;
    return map_matrix_status(numelyx_core_matrix_inverse(matrix, inverse));
}

numelyx_status numelyx_matrix_solve(
    const numelyx_matrix *matrix,
    const numelyx_matrix *right_hand_side,
    numelyx_matrix **solution
)
{
    if (solution == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *solution = NULL;
    return map_matrix_status(numelyx_core_matrix_solve(
        matrix,
        right_hand_side,
        solution
    ));
}



numelyx_status numelyx_matrix_get_storage_kind(
    const numelyx_matrix *matrix,
    numelyx_storage_kind *kind
)
{
    numelyx_storage_kind mapped;

    if (matrix == NULL || kind == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }

    switch (numelyx_core_matrix_storage_kind(matrix)) {
        case NUMERUS_STORAGE_DENSE: mapped = NUMELYX_STORAGE_DENSE; break;
        case NUMERUS_STORAGE_UPPER_TRIANGULAR: mapped = NUMELYX_STORAGE_UPPER_TRIANGULAR; break;
        case NUMERUS_STORAGE_LOWER_TRIANGULAR: mapped = NUMELYX_STORAGE_LOWER_TRIANGULAR; break;
        case NUMERUS_STORAGE_DIAGONAL: mapped = NUMELYX_STORAGE_DIAGONAL; break;
        case NUMERUS_STORAGE_IDENTITY: mapped = NUMELYX_STORAGE_IDENTITY; break;
        case NUMERUS_STORAGE_CONSTANT: mapped = NUMELYX_STORAGE_CONSTANT; break;
        case NUMERUS_STORAGE_ZERO: mapped = NUMELYX_STORAGE_ZERO; break;
        case NUMERUS_STORAGE_SCALED_IDENTITY: mapped = NUMELYX_STORAGE_SCALED_IDENTITY; break;
        case NUMERUS_STORAGE_SPARSE: mapped = NUMELYX_STORAGE_SPARSE; break;
        case NUMERUS_STORAGE_SYMMETRIC: mapped = NUMELYX_STORAGE_SYMMETRIC; break;
        case NUMERUS_STORAGE_BANDED: mapped = NUMELYX_STORAGE_BANDED; break;
        default: return NUMELYX_NUMERICAL_ERROR;
    }

    *kind = mapped;
    return NUMELYX_SUCCESS;
}

numelyx_status numelyx_matrix_get_flags(
    numelyx_matrix *matrix,
    numelyx_matrix_flags *flags
)
{
    numelyx_core_matrix_flags core_flags;
    numelyx_matrix_flags result;
    numelyx_core_matrix_status status;

    if (flags == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }

    status = numelyx_core_matrix_get_flags(matrix, &core_flags);
    if (status != NUMERUS_MATRIX_SUCCESS) {
        return map_matrix_status(status);
    }

    result.square = core_flags.square ? 1 : 0;
    result.zero = core_flags.zero ? 1 : 0;
    result.diagonal = core_flags.diagonal ? 1 : 0;
    result.upper_triangular = core_flags.upper_triangular ? 1 : 0;
    result.lower_triangular = core_flags.lower_triangular ? 1 : 0;
    result.symmetric = core_flags.symmetric ? 1 : 0;
    result.identity = core_flags.identity ? 1 : 0;
    *flags = result;
    return NUMELYX_SUCCESS;
}

numelyx_status numelyx_matrix_log_determinant(
    const numelyx_matrix *matrix,
    int *sign,
    double *log_abs_determinant
)
{
    return map_matrix_status(numelyx_core_matrix_log_determinant(
        matrix, sign, log_abs_determinant));
}

numelyx_status numelyx_matrix_qr_decompose(
    const numelyx_matrix *matrix,
    numelyx_matrix **q,
    numelyx_matrix **r
)
{
    if (q == NULL || r == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *q = NULL;
    *r = NULL;
    if (q == r) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    return map_matrix_status(numelyx_core_matrix_qr_decompose(matrix, q, r));
}

numelyx_status numelyx_matrix_svd(
    const numelyx_matrix *matrix,
    numelyx_matrix **u,
    numelyx_matrix **singular_values,
    numelyx_matrix **vt
)
{
    if (u == NULL || singular_values == NULL || vt == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *u = NULL;
    *singular_values = NULL;
    *vt = NULL;
    if (u == singular_values || u == vt || singular_values == vt) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    return map_matrix_status(numelyx_core_matrix_svd(
        matrix, u, singular_values, vt));
}

numelyx_status numelyx_matrix_pseudoinverse(
    const numelyx_matrix *matrix,
    numelyx_matrix **pseudoinverse
)
{
    if (pseudoinverse == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *pseudoinverse = NULL;
    return map_matrix_status(numelyx_core_matrix_pseudoinverse(
        matrix, pseudoinverse));
}

numelyx_status numelyx_matrix_cholesky(
    const numelyx_matrix *matrix,
    numelyx_matrix **lower
)
{
    if (lower == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *lower = NULL;
    return map_matrix_status(numelyx_core_matrix_cholesky(matrix, lower));
}

numelyx_status numelyx_matrix_least_squares(
    const numelyx_matrix *matrix,
    const numelyx_matrix *right_hand_side,
    numelyx_matrix **solution
)
{
    if (solution == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *solution = NULL;
    return map_matrix_status(numelyx_core_matrix_least_squares(
        matrix, right_hand_side, solution));
}

numelyx_status numelyx_matrix_weighted_least_squares(
    const numelyx_matrix *matrix,
    const numelyx_matrix *right_hand_side,
    const double *weights,
    size_t weights_count,
    numelyx_matrix **solution
)
{
    if (solution == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *solution = NULL;
    return map_matrix_status(numelyx_core_matrix_weighted_least_squares(
        matrix, right_hand_side, weights, weights_count, solution));
}

numelyx_status numelyx_matrix_rank(
    const numelyx_matrix *matrix,
    size_t *rank
)
{
    return map_matrix_status(numelyx_core_matrix_rank(matrix, rank));
}

numelyx_status numelyx_matrix_solve_triangular(
    const numelyx_matrix *triangular,
    const numelyx_matrix *right_hand_side,
    int lower,
    int transpose,
    numelyx_matrix **solution
)
{
    if (solution == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *solution = NULL;
    return map_matrix_status(numelyx_core_matrix_solve_triangular(
        triangular, right_hand_side, lower != 0, transpose != 0, solution));
}

numelyx_status numelyx_matrix_norm_frobenius(
    const numelyx_matrix *matrix,
    double *norm
)
{
    return map_matrix_status(numelyx_core_matrix_norm_frobenius(matrix, norm));
}

numelyx_status numelyx_matrix_norm_one(
    const numelyx_matrix *matrix,
    double *norm
)
{
    return map_matrix_status(numelyx_core_matrix_norm_one(matrix, norm));
}

numelyx_status numelyx_matrix_norm_infinity(
    const numelyx_matrix *matrix,
    double *norm
)
{
    return map_matrix_status(numelyx_core_matrix_norm_infinity(matrix, norm));
}

numelyx_status numelyx_matrix_trace(
    const numelyx_matrix *matrix,
    double *trace
)
{
    return map_matrix_status(numelyx_core_matrix_trace(matrix, trace));
}


numelyx_status numelyx_matrix_symmetric_eigen(
    const numelyx_matrix *matrix,
    numelyx_matrix **eigenvalues,
    numelyx_matrix **eigenvectors
)
{
    if (eigenvalues == NULL || eigenvectors == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *eigenvalues = NULL;
    *eigenvectors = NULL;
    if (eigenvalues == eigenvectors) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    return map_matrix_status(numelyx_core_matrix_symmetric_eigen(
        matrix, eigenvalues, eigenvectors));
}

numelyx_status numelyx_matrix_ldlt_decompose(
    const numelyx_matrix *matrix,
    numelyx_matrix **lower,
    numelyx_matrix **diagonal
)
{
    if (lower == NULL || diagonal == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *lower = NULL;
    *diagonal = NULL;
    if (lower == diagonal) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    return map_matrix_status(numelyx_core_matrix_ldlt_decompose(
        matrix, lower, diagonal));
}

numelyx_status numelyx_matrix_spectral_norm(
    const numelyx_matrix *matrix,
    double *norm
)
{
    return map_matrix_status(numelyx_core_matrix_spectral_norm(matrix, norm));
}

numelyx_status numelyx_matrix_condition_estimate_one(
    const numelyx_matrix *matrix,
    double *condition_estimate
)
{
    return map_matrix_status(numelyx_core_matrix_condition_estimate_one(
        matrix, condition_estimate));
}

numelyx_status numelyx_matrix_classify_definiteness(
    const numelyx_matrix *matrix,
    int *positive_definite,
    int *positive_semidefinite
)
{
    return map_matrix_status(numelyx_core_matrix_classify_definiteness(
        matrix, positive_definite, positive_semidefinite));
}

numelyx_status numelyx_matrix_classify_system(
    const numelyx_matrix *matrix,
    const numelyx_matrix *right_hand_side,
    numelyx_matrix_solution_kind *classification
)
{
    numelyx_core_matrix_solution_kind core_classification;
    numelyx_core_matrix_status status;

    if (classification == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    status = numelyx_core_matrix_classify_system(
        matrix, right_hand_side, &core_classification);
    if (status != NUMERUS_MATRIX_SUCCESS) {
        return map_matrix_status(status);
    }

    switch (core_classification) {
        case NUMERUS_MATRIX_SOLUTION_UNIQUE:
            *classification = NUMELYX_MATRIX_SOLUTION_UNIQUE;
            break;
        case NUMERUS_MATRIX_SOLUTION_INFINITE:
            *classification = NUMELYX_MATRIX_SOLUTION_INFINITE;
            break;
        case NUMERUS_MATRIX_SOLUTION_INCONSISTENT:
            *classification = NUMELYX_MATRIX_SOLUTION_INCONSISTENT;
            break;
        default:
            return NUMELYX_NUMERICAL_ERROR;
    }
    return NUMELYX_SUCCESS;
}

numelyx_status numelyx_matrix_sum(
    const numelyx_matrix *matrix,
    double *sum
)
{
    return map_matrix_status(numelyx_core_matrix_sum(matrix, sum));
}

numelyx_status numelyx_matrix_mean(
    const numelyx_matrix *matrix,
    double *mean
)
{
    return map_matrix_status(numelyx_core_matrix_mean(matrix, mean));
}

numelyx_status numelyx_matrix_variance(
    const numelyx_matrix *matrix,
    int sample,
    double *variance
)
{
    if (sample != 0 && sample != 1) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    return map_matrix_status(numelyx_core_matrix_variance(
        matrix, sample != 0, variance));
}

numelyx_status numelyx_matrix_standard_deviation(
    const numelyx_matrix *matrix,
    int sample,
    double *standard_deviation
)
{
    if (sample != 0 && sample != 1) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    return map_matrix_status(numelyx_core_matrix_standard_deviation(
        matrix, sample != 0, standard_deviation));
}

numelyx_status numelyx_matrix_min(
    const numelyx_matrix *matrix,
    double *minimum
)
{
    return map_matrix_status(numelyx_core_matrix_min(matrix, minimum));
}

numelyx_status numelyx_matrix_max(
    const numelyx_matrix *matrix,
    double *maximum
)
{
    return map_matrix_status(numelyx_core_matrix_max(matrix, maximum));
}

numelyx_status numelyx_matrix_retain(numelyx_matrix *matrix)
{
    return map_matrix_status(numelyx_core_matrix_retain(matrix));
}

void numelyx_matrix_release(numelyx_matrix *matrix)
{
    numelyx_core_matrix_destroy(matrix);
}

size_t numelyx_matrix_rows(const numelyx_matrix *matrix)
{
    return numelyx_core_matrix_rows(matrix);
}

size_t numelyx_matrix_columns(const numelyx_matrix *matrix)
{
    return numelyx_core_matrix_columns(matrix);
}

numelyx_status numelyx_matrix_get(
    const numelyx_matrix *matrix,
    size_t row,
    size_t column,
    double *value
)
{
    return map_matrix_status(numelyx_core_matrix_get(
        matrix,
        row,
        column,
        value
    ));
}
