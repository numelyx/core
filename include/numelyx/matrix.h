#ifndef NUMELYX_MATRIX_H
#define NUMELYX_MATRIX_H

#include <stddef.h>

#include <numelyx/status.h>

NUMELYX_BEGIN_DECLS

/** Opaque Matrix handle. The implementation and storage layout are private. */
typedef struct numelyx_core_matrix numelyx_matrix;

/**
 * Create a dense row-major Matrix. The input values are copied.
 *
 * rows and columns must be positive and values must point to at least
 * rows * columns doubles. On failure, *matrix is set to NULL.
 */
NUMELYX_API numelyx_status numelyx_matrix_create_dense(
    size_t rows,
    size_t columns,
    const double *values,
    numelyx_matrix **matrix
);

/** Create an identity Matrix without allocating an element buffer. */
NUMELYX_API numelyx_status numelyx_matrix_create_identity(
    size_t size,
    numelyx_matrix **matrix
);

/** Create an all-zero Matrix without allocating an element buffer. */
NUMELYX_API numelyx_status numelyx_matrix_create_zero(
    size_t rows,
    size_t columns,
    numelyx_matrix **matrix
);

/** Create a constant-valued Matrix using one stored scalar. */
NUMELYX_API numelyx_status numelyx_matrix_create_constant(
    size_t rows,
    size_t columns,
    double value,
    numelyx_matrix **matrix
);

/** Create a diagonal Matrix from a copied diagonal-value buffer. */
NUMELYX_API numelyx_status numelyx_matrix_create_diagonal(
    size_t size,
    const double *diagonal,
    numelyx_matrix **matrix
);

/** Retain one additional owned reference; the caller must already own one. */
NUMELYX_API numelyx_status numelyx_matrix_retain(numelyx_matrix *matrix);

/** Release one owned reference. NULL is accepted. */
NUMELYX_API void numelyx_matrix_release(numelyx_matrix *matrix);

/** Return dimensions, or zero when matrix is NULL. */
NUMELYX_API size_t numelyx_matrix_rows(const numelyx_matrix *matrix);
NUMELYX_API size_t numelyx_matrix_columns(const numelyx_matrix *matrix);

/** Read an element. The output value is unchanged on failure. */
NUMELYX_API numelyx_status numelyx_matrix_get(
    const numelyx_matrix *matrix,
    size_t row,
    size_t column,
    double *value
);

NUMELYX_END_DECLS

#endif /* NUMELYX_MATRIX_H */
