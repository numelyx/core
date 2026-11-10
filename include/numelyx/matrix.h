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

/** Create a scaled identity Matrix with one stored diagonal scalar. */
NUMELYX_API numelyx_status numelyx_matrix_create_scaled_identity(
    size_t size,
    double diagonal_value,
    numelyx_matrix **matrix
);

/** Create packed upper-triangular storage (row-wise, on and above diagonal). */
NUMELYX_API numelyx_status numelyx_matrix_create_upper_triangular(
    size_t size,
    const double *packed_values,
    numelyx_matrix **matrix
);

/** Create packed lower-triangular storage (row-wise, on and below diagonal). */
NUMELYX_API numelyx_status numelyx_matrix_create_lower_triangular(
    size_t size,
    const double *packed_values,
    numelyx_matrix **matrix
);

/** Create symmetric storage from packed lower-triangular values. */
/** One sparse override, addressed by row-major linear index. */
typedef struct {
    size_t index;
    double value;
} numelyx_sparse_entry;

/** Create sparse storage with a default value and unique linear-index overrides. */
NUMELYX_API numelyx_status numelyx_matrix_create_sparse(
    size_t rows,
    size_t columns,
    double default_value,
    const numelyx_sparse_entry *entries,
    size_t count,
    numelyx_matrix **matrix
);

/** Create fixed-bandwidth storage from row-major packed band rows.
 * Each row stores lower_bandwidth + upper_bandwidth + 1 values, with the
 * main diagonal at slot upper_bandwidth; out-of-matrix slots are ignored.
 */
NUMELYX_API numelyx_status numelyx_matrix_create_banded(
    size_t rows,
    size_t columns,
    size_t lower_bandwidth,
    size_t upper_bandwidth,
    const double *packed_band_values,
    numelyx_matrix **matrix
);


NUMELYX_API numelyx_status numelyx_matrix_create_symmetric(
    size_t size,
    const double *packed_lower_values,
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

/** Create a lazy view with rows in reverse order. */
NUMELYX_API numelyx_status numelyx_matrix_flip_rows(
    numelyx_matrix *parent,
    numelyx_matrix **view
);

/** Create a lazy view with columns in reverse order. */
NUMELYX_API numelyx_status numelyx_matrix_flip_columns(
    numelyx_matrix *parent,
    numelyx_matrix **view
);

/** Select rows by copied zero-based indices; repeated indices are allowed. */
NUMELYX_API numelyx_status numelyx_matrix_select_rows(
    numelyx_matrix *parent,
    const size_t *indices,
    size_t count,
    numelyx_matrix **view
);

/** Select columns by copied zero-based indices; repeated indices are allowed. */
NUMELYX_API numelyx_status numelyx_matrix_select_columns(
    numelyx_matrix *parent,
    const size_t *indices,
    size_t count,
    numelyx_matrix **view
);

/** Create a lazy row-major reshape with the same element count. */
NUMELYX_API numelyx_status numelyx_matrix_reshape(
    numelyx_matrix *parent,
    size_t rows,
    size_t columns,
    numelyx_matrix **view
);

/** Create a lazy 1 × N row-vector view of all parent values in row-major order. */
NUMELYX_API numelyx_status numelyx_matrix_flatten(
    numelyx_matrix *parent,
    numelyx_matrix **view
);

/** Remove one row or column as a lazy view. */
NUMELYX_API numelyx_status numelyx_matrix_remove_row(
    numelyx_matrix *parent, size_t row, numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_remove_column(
    numelyx_matrix *parent, size_t column, numelyx_matrix **view
);

/** Swap two rows or columns as a lazy view. */
NUMELYX_API numelyx_status numelyx_matrix_swap_rows(
    numelyx_matrix *parent, size_t row1, size_t row2, numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_swap_columns(
    numelyx_matrix *parent, size_t column1, size_t column2, numelyx_matrix **view
);

/** Create lazy 90-degree clockwise, 180-degree, or counter-clockwise rotations. */
NUMELYX_API numelyx_status numelyx_matrix_rotate_90_clockwise(
    numelyx_matrix *parent, numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_rotate_180(
    numelyx_matrix *parent, numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_rotate_90_counterclockwise(
    numelyx_matrix *parent, numelyx_matrix **view
);

/** Join two matrices horizontally or vertically without copying values. */
NUMELYX_API numelyx_status numelyx_matrix_join_horizontal(
    numelyx_matrix *left, numelyx_matrix *right, numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_join_vertical(
    numelyx_matrix *top, numelyx_matrix *bottom, numelyx_matrix **view
);

/** Add constant or zero padding as a lazy view. */
NUMELYX_API numelyx_status numelyx_matrix_pad(
    numelyx_matrix *parent,
    size_t top, size_t bottom, size_t left, size_t right,
    double value, numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_zero_extend(
    numelyx_matrix *parent,
    size_t top, size_t bottom, size_t left, size_t right,
    numelyx_matrix **view
);

/** Repeat a matrix lazily by row and column factors. */
NUMELYX_API numelyx_status numelyx_matrix_repeat(
    numelyx_matrix *parent,
    size_t row_repetitions,
    size_t column_repetitions,
    numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_block_diagonal(
    numelyx_matrix *parent,
    size_t repetitions,
    numelyx_matrix **view
);

/** Assemble a row-major block grid; the pointer array is borrowed during call. */
NUMELYX_API numelyx_status numelyx_matrix_block_grid(
    numelyx_matrix *const *blocks,
    size_t block_row_count,
    size_t block_column_count,
    numelyx_matrix **view
);

/** Extract a diagonal vector or create a square diagonal view from a vector. */
NUMELYX_API numelyx_status numelyx_matrix_diagonal_extract(
    numelyx_matrix *parent,
    ptrdiff_t offset,
    numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_diagonal_from_vector(
    numelyx_matrix *vector,
    ptrdiff_t offset,
    numelyx_matrix **view
);

/** Apply full row/column permutations; index arrays are copied. */
NUMELYX_API numelyx_status numelyx_matrix_permute_rows(
    numelyx_matrix *parent,
    const size_t *permutation,
    size_t count,
    numelyx_matrix **view
);
NUMELYX_API numelyx_status numelyx_matrix_permute_columns(
    numelyx_matrix *parent,
    const size_t *permutation,
    size_t count,
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

/** Compute the determinant of a square Matrix. */
NUMELYX_API numelyx_status numelyx_matrix_determinant(
    numelyx_matrix *matrix,
    double *determinant
);

/** Compute an independent dense inverse using the Core LU-based path. */
NUMELYX_API numelyx_status numelyx_matrix_inverse(
    const numelyx_matrix *matrix,
    numelyx_matrix **inverse
);

/** Solve A X = B using one LU factorization for all right-hand sides. */
NUMELYX_API numelyx_status numelyx_matrix_solve(
    const numelyx_matrix *matrix,
    const numelyx_matrix *right_hand_side,
    numelyx_matrix **solution
);

/** Stable public identifiers for the root storage representation. */
typedef enum {
    NUMELYX_STORAGE_DENSE = 0,
    NUMELYX_STORAGE_UPPER_TRIANGULAR,
    NUMELYX_STORAGE_LOWER_TRIANGULAR,
    NUMELYX_STORAGE_DIAGONAL,
    NUMELYX_STORAGE_IDENTITY,
    NUMELYX_STORAGE_CONSTANT,
    NUMELYX_STORAGE_ZERO,
    NUMELYX_STORAGE_SCALED_IDENTITY,
    NUMELYX_STORAGE_SPARSE,
    NUMELYX_STORAGE_SYMMETRIC,
    NUMELYX_STORAGE_BANDED
} numelyx_storage_kind;

/** Logical structural properties, computed/cached by the Core Matrix layer. */
typedef struct {
    int square;
    int zero;
    int diagonal;
    int upper_triangular;
    int lower_triangular;
    int symmetric;
    int identity;
} numelyx_matrix_flags;

/** Return the underlying root storage kind; views report their root kind. */
NUMELYX_API numelyx_status numelyx_matrix_get_storage_kind(
    const numelyx_matrix *matrix, numelyx_storage_kind *kind
);

/** Compute or retrieve logical structural flags; output is unchanged on error. */
NUMELYX_API numelyx_status numelyx_matrix_get_flags(
    numelyx_matrix *matrix, numelyx_matrix_flags *flags
);

/** Compute signed log(|det(A)|) without forming the determinant product. */
NUMELYX_API numelyx_status numelyx_matrix_log_determinant(
    const numelyx_matrix *matrix, int *sign, double *log_abs_determinant
);

/** Reduced Householder QR decomposition A = Q R. */
NUMELYX_API numelyx_status numelyx_matrix_qr_decompose(
    const numelyx_matrix *matrix, numelyx_matrix **q, numelyx_matrix **r
);

/** Reduced SVD A = U Sigma V^T; singular values are descending. */
NUMELYX_API numelyx_status numelyx_matrix_svd(
    const numelyx_matrix *matrix,
    numelyx_matrix **u,
    numelyx_matrix **singular_values,
    numelyx_matrix **vt
);

/** Moore-Penrose pseudoinverse computed from the reduced SVD. */
NUMELYX_API numelyx_status numelyx_matrix_pseudoinverse(
    const numelyx_matrix *matrix, numelyx_matrix **pseudoinverse
);

/** Cholesky factor L for symmetric positive-definite A = L L^T. */
NUMELYX_API numelyx_status numelyx_matrix_cholesky(
    const numelyx_matrix *matrix, numelyx_matrix **lower
);

/** Solve min ||A X - B||_2 using column-pivoted QR. */
NUMELYX_API numelyx_status numelyx_matrix_least_squares(
    const numelyx_matrix *matrix,
    const numelyx_matrix *right_hand_side,
    numelyx_matrix **solution
);

/** Solve weighted least squares using nonnegative row weights. */
NUMELYX_API numelyx_status numelyx_matrix_weighted_least_squares(
    const numelyx_matrix *matrix,
    const numelyx_matrix *right_hand_side,
    const double *weights,
    size_t weights_count,
    numelyx_matrix **solution
);

/** Estimate numerical rank using scale-aware complete pivoting. */
NUMELYX_API numelyx_status numelyx_matrix_rank(
    const numelyx_matrix *matrix, size_t *rank
);

/** Solve T X = B or T^T X = B by triangular substitution. */
NUMELYX_API numelyx_status numelyx_matrix_solve_triangular(
    const numelyx_matrix *triangular,
    const numelyx_matrix *right_hand_side,
    int lower,
    int transpose,
    numelyx_matrix **solution
);

/** Compute the Frobenius, induced 1-, or induced infinity-norm. */
NUMELYX_API numelyx_status numelyx_matrix_norm_frobenius(
    const numelyx_matrix *matrix, double *norm
);
NUMELYX_API numelyx_status numelyx_matrix_norm_one(
    const numelyx_matrix *matrix, double *norm
);
NUMELYX_API numelyx_status numelyx_matrix_norm_infinity(
    const numelyx_matrix *matrix, double *norm
);

/** Compute the trace of a square matrix. */
NUMELYX_API numelyx_status numelyx_matrix_trace(
    const numelyx_matrix *matrix, double *trace
);

/** Numerical classification of solutions to A X = B. */
typedef enum {
    NUMELYX_MATRIX_SOLUTION_UNIQUE = 0,
    NUMELYX_MATRIX_SOLUTION_INFINITE,
    NUMELYX_MATRIX_SOLUTION_INCONSISTENT
} numelyx_matrix_solution_kind;

/** Symmetric eigendecomposition; eigenvalues are descending, eigenvectors are columns. */
NUMELYX_API numelyx_status numelyx_matrix_symmetric_eigen(
    const numelyx_matrix *matrix,
    numelyx_matrix **eigenvalues,
    numelyx_matrix **eigenvectors
);

/** Unpivoted LDL^T factorization of a symmetric matrix. */
NUMELYX_API numelyx_status numelyx_matrix_ldlt_decompose(
    const numelyx_matrix *matrix,
    numelyx_matrix **lower,
    numelyx_matrix **diagonal
);

/** Spectral norm and 1-norm condition estimate. */
NUMELYX_API numelyx_status numelyx_matrix_spectral_norm(
    const numelyx_matrix *matrix, double *norm
);
NUMELYX_API numelyx_status numelyx_matrix_condition_estimate_one(
    const numelyx_matrix *matrix, double *condition_estimate
);

/** Classify positive definiteness and semidefiniteness; outputs must differ. */
NUMELYX_API numelyx_status numelyx_matrix_classify_definiteness(
    const numelyx_matrix *matrix,
    int *positive_definite,
    int *positive_semidefinite
);

/** Classify the solution set of A X = B, including multiple RHS columns. */
NUMELYX_API numelyx_status numelyx_matrix_classify_system(
    const numelyx_matrix *matrix,
    const numelyx_matrix *right_hand_side,
    numelyx_matrix_solution_kind *classification
);

/** Basic reductions. Variance uses sample=0 for population and sample=1 for sample variance. */
NUMELYX_API numelyx_status numelyx_matrix_sum(
    const numelyx_matrix *matrix, double *sum
);
NUMELYX_API numelyx_status numelyx_matrix_mean(
    const numelyx_matrix *matrix, double *mean
);
NUMELYX_API numelyx_status numelyx_matrix_variance(
    const numelyx_matrix *matrix, int sample, double *variance
);
NUMELYX_API numelyx_status numelyx_matrix_standard_deviation(
    const numelyx_matrix *matrix, int sample, double *standard_deviation
);
NUMELYX_API numelyx_status numelyx_matrix_min(
    const numelyx_matrix *matrix, double *minimum
);
NUMELYX_API numelyx_status numelyx_matrix_max(
    const numelyx_matrix *matrix, double *maximum
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
