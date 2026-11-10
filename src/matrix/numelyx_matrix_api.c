#include <numelyx/matrix.h>

#include "numelyx_core_matrix_core.h"
#include "numelyx_core_matrix_types.h"
#include "numelyx_core_matrix_views.h"

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
