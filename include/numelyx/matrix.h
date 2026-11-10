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

/**
 * Create a lazy transpose view. The view retains its parent and does not copy
 * matrix values. Release the parent reference independently after creation.
 */
NUMELYX_API numelyx_status numelyx_matrix_transpose(
    numelyx_matrix *parent,
    numelyx_matrix **view
);

/**
 * Create a lazy rectangular slice. Counts must be positive and the selected
 * range must fit inside the parent. The view retains its parent.
 */
NUMELYX_API numelyx_status numelyx_matrix_slice(
    numelyx_matrix *parent,
    size_t row_start,
    size_t row_count,
    size_t column_start,
    size_t column_count,
    numelyx_matrix **view
);

/** Create a lazy element-wise sum; the result retains both parents. */
NUMELYX_API numelyx_status numelyx_matrix_add(
    numelyx_matrix *left,
    numelyx_matrix *right,
    numelyx_matrix **result
);

/** Create a lazy element-wise difference; the result retains both parents. */
NUMELYX_API numelyx_status numelyx_matrix_subtract(
    numelyx_matrix *left,
    numelyx_matrix *right,
    numelyx_matrix **result
);

/** Create a lazy Hadamard (element-wise) product. */
NUMELYX_API numelyx_status numelyx_matrix_hadamard_product(
    numelyx_matrix *left,
    numelyx_matrix *right,
    numelyx_matrix **result
);

/** Multiply left × right into an independent dense Matrix. */
NUMELYX_API numelyx_status numelyx_matrix_multiply(
    const numelyx_matrix *left,
    const numelyx_matrix *right,
    numelyx_matrix **result
);

/** Create a lazy scalar multiplication view. */
NUMELYX_API numelyx_status numelyx_matrix_scale(
    numelyx_matrix *source,
    double scalar,
    numelyx_matrix **result
);

/** Materialize any Matrix/view into independent dense storage. */
NUMELYX_API numelyx_status numelyx_matrix_materialize(
    const numelyx_matrix *source,
    numelyx_matrix **result
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
