#include <numelyx/matrix.h>

#include "numelyx_core_matrix_core.h"
#include "numelyx_core_matrix_types.h"

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

numelyx_status numelyx_matrix_create_dense(
    size_t rows,
    size_t columns,
    const double *values,
    numelyx_matrix **matrix
)
{
    numelyx_core_matrix_data data = {0};
    int status;

    if (matrix == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }

    *matrix = NULL;
    if (values == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }

    data.values = values;
    status = numelyx_core_matrix_create(
        NUMERUS_STORAGE_DENSE,
        rows,
        columns,
        &data,
        matrix
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
