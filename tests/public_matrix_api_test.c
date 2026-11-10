#include <math.h>

#include <numelyx/matrix.h>

int main(void)
{
    double values[] = {1.0, 2.0, 3.0, 4.0};
    double value = -1.0;
    numelyx_matrix *matrix = NULL;

    if (numelyx_matrix_create_dense(2, 2, values, &matrix) != NUMELYX_SUCCESS) {
        return 1;
    }
    values[2] = 99.0;
    if (matrix == NULL
        || numelyx_matrix_rows(matrix) != 2
        || numelyx_matrix_columns(matrix) != 2) {
        numelyx_matrix_release(matrix);
        return 2;
    }
    if (numelyx_matrix_get(matrix, 1, 0, &value) != NUMELYX_SUCCESS
        || value != 3.0) {
        numelyx_matrix_release(matrix);
        return 3;
    }

    value = 17.0;
    if (numelyx_matrix_get(matrix, 2, 0, &value) != NUMELYX_OUT_OF_BOUNDS
        || value != 17.0) {
        numelyx_matrix_release(matrix);
        return 4;
    }
    if (numelyx_matrix_retain(matrix) != NUMELYX_SUCCESS) {
        numelyx_matrix_release(matrix);
        return 5;
    }

    numelyx_matrix_release(matrix);
    numelyx_matrix_release(matrix);

    if (numelyx_matrix_create_dense(0, 2, values, &matrix)
            != NUMELYX_INVALID_ARGUMENT
        || matrix != NULL) {
        numelyx_matrix_release(matrix);
        return 6;
    }
    if (numelyx_matrix_create_dense(2, 2, NULL, &matrix)
            != NUMELYX_INVALID_ARGUMENT
        || matrix != NULL) {
        numelyx_matrix_release(matrix);
        return 7;
    }
    if (numelyx_matrix_retain(NULL) != NUMELYX_INVALID_ARGUMENT) {
        return 8;
    }

    {
        const double source_values[] = {1.0, 2.0, 3.0, 4.0};
        numelyx_matrix *parent = NULL;
        numelyx_matrix *transpose = NULL;
        numelyx_matrix *slice = NULL;

        if (numelyx_matrix_create_dense(2, 2, source_values, &parent)
                != NUMELYX_SUCCESS
            || numelyx_matrix_transpose(parent, &transpose) != NUMELYX_SUCCESS
            || numelyx_matrix_slice(transpose, 0, 1, 1, 1, &slice)
                != NUMELYX_SUCCESS) {
            numelyx_matrix_release(slice);
            numelyx_matrix_release(transpose);
            numelyx_matrix_release(parent);
            return 13;
        }

        numelyx_matrix_release(parent);
        numelyx_matrix_release(transpose);
        if (numelyx_matrix_rows(slice) != 1
            || numelyx_matrix_columns(slice) != 1
            || numelyx_matrix_get(slice, 0, 0, &value) != NUMELYX_SUCCESS
            || value != 3.0) {
            numelyx_matrix_release(slice);
            return 14;
        }
        numelyx_matrix_release(slice);
    }

    {
        const double left_values[] = {1.0, 2.0, 3.0, 4.0};
        const double right_values[] = {5.0, 6.0, 7.0, 8.0};
        numelyx_matrix *left = NULL;
        numelyx_matrix *right = NULL;
        numelyx_matrix *sum = NULL;
        numelyx_matrix *product = NULL;
        numelyx_matrix *hadamard = NULL;
        numelyx_matrix *scaled = NULL;
        numelyx_matrix *materialized = NULL;

        if (numelyx_matrix_create_dense(2, 2, left_values, &left)
                != NUMELYX_SUCCESS
            || numelyx_matrix_create_dense(2, 2, right_values, &right)
                != NUMELYX_SUCCESS
            || numelyx_matrix_add(left, right, &sum) != NUMELYX_SUCCESS
            || numelyx_matrix_multiply(left, right, &product) != NUMELYX_SUCCESS
            || numelyx_matrix_hadamard_product(left, right, &hadamard)
                != NUMELYX_SUCCESS
            || numelyx_matrix_scale(left, -2.0, &scaled) != NUMELYX_SUCCESS
            || numelyx_matrix_materialize(sum, &materialized) != NUMELYX_SUCCESS) {
            numelyx_matrix_release(materialized);
            numelyx_matrix_release(scaled);
            numelyx_matrix_release(hadamard);
            numelyx_matrix_release(product);
            numelyx_matrix_release(sum);
            numelyx_matrix_release(right);
            numelyx_matrix_release(left);
            return 15;
        }

        /* Lazy results own their parents; eager multiplication does not. */
        numelyx_matrix_release(left);
        left = NULL;
        numelyx_matrix_release(right);
        right = NULL;

        if (numelyx_matrix_get(sum, 1, 1, &value) != NUMELYX_SUCCESS
                || value != 12.0
            || numelyx_matrix_get(product, 1, 0, &value) != NUMELYX_SUCCESS
                || value != 43.0
            || numelyx_matrix_get(hadamard, 1, 1, &value) != NUMELYX_SUCCESS
                || value != 32.0
            || numelyx_matrix_get(scaled, 0, 1, &value) != NUMELYX_SUCCESS
                || value != -4.0
            || numelyx_matrix_get(materialized, 0, 1, &value) != NUMELYX_SUCCESS
                || value != 8.0) {
            numelyx_matrix_release(materialized);
            numelyx_matrix_release(scaled);
            numelyx_matrix_release(hadamard);
            numelyx_matrix_release(product);
            numelyx_matrix_release(sum);
            numelyx_matrix_release(right);
            numelyx_matrix_release(left);
            return 16;
        }

        numelyx_matrix_release(materialized);
        numelyx_matrix_release(scaled);
        numelyx_matrix_release(hadamard);
        numelyx_matrix_release(product);
        numelyx_matrix_release(sum);
        numelyx_matrix_release(right);
        numelyx_matrix_release(left);
    }

    {
        const double source_values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        const size_t row_indices[] = {1, 0};
        const size_t column_indices[] = {2, 0};
        numelyx_matrix *parent = NULL;
        numelyx_matrix *rows = NULL;
        numelyx_matrix *columns = NULL;
        numelyx_matrix *reshaped = NULL;
        numelyx_matrix *flattened = NULL;
        numelyx_matrix *flipped_rows = NULL;
        numelyx_matrix *flipped_columns = NULL;

        if (numelyx_matrix_create_dense(2, 3, source_values, &parent)
                != NUMELYX_SUCCESS
            || numelyx_matrix_select_rows(parent, row_indices, 2, &rows)
                != NUMELYX_SUCCESS
            || numelyx_matrix_select_columns(parent, column_indices, 2, &columns)
                != NUMELYX_SUCCESS
            || numelyx_matrix_reshape(parent, 3, 2, &reshaped)
                != NUMELYX_SUCCESS
            || numelyx_matrix_flatten(parent, &flattened) != NUMELYX_SUCCESS
            || numelyx_matrix_flip_rows(parent, &flipped_rows) != NUMELYX_SUCCESS
            || numelyx_matrix_flip_columns(parent, &flipped_columns)
                != NUMELYX_SUCCESS) {
            numelyx_matrix_release(flipped_columns);
            numelyx_matrix_release(flipped_rows);
            numelyx_matrix_release(flattened);
            numelyx_matrix_release(reshaped);
            numelyx_matrix_release(columns);
            numelyx_matrix_release(rows);
            numelyx_matrix_release(parent);
            return 19;
        }

        numelyx_matrix_release(parent);
        if (numelyx_matrix_get(rows, 0, 2, &value) != NUMELYX_SUCCESS
                || value != 6.0
            || numelyx_matrix_get(columns, 1, 1, &value) != NUMELYX_SUCCESS
                || value != 4.0
            || numelyx_matrix_get(reshaped, 1, 0, &value) != NUMELYX_SUCCESS
                || value != 3.0
            || numelyx_matrix_get(flattened, 0, 4, &value) != NUMELYX_SUCCESS
                || value != 5.0
            || numelyx_matrix_get(flipped_rows, 0, 0, &value) != NUMELYX_SUCCESS
                || value != 4.0
            || numelyx_matrix_get(flipped_columns, 0, 0, &value)
                != NUMELYX_SUCCESS
            || value != 3.0) {
            numelyx_matrix_release(flipped_columns);
            numelyx_matrix_release(flipped_rows);
            numelyx_matrix_release(flattened);
            numelyx_matrix_release(reshaped);
            numelyx_matrix_release(columns);
            numelyx_matrix_release(rows);
            return 20;
        }

        numelyx_matrix_release(flipped_columns);
        numelyx_matrix_release(flipped_rows);
        numelyx_matrix_release(flattened);
        numelyx_matrix_release(reshaped);
        numelyx_matrix_release(columns);
        numelyx_matrix_release(rows);
    }

    {
        const double source_values[] = {1.0, 2.0, 3.0, 4.0};
        const double vector_values[] = {5.0, 6.0};
        const size_t permutation[] = {1, 0};
        numelyx_matrix *parent = NULL;
        numelyx_matrix *rotated = NULL;
        numelyx_matrix *rotated_180 = NULL;
        numelyx_matrix *rotated_ccw = NULL;
        numelyx_matrix *removed = NULL;
        numelyx_matrix *swapped = NULL;
        numelyx_matrix *swapped_rows = NULL;
        numelyx_matrix *permuted_columns = NULL;
        numelyx_matrix *joined = NULL;
        numelyx_matrix *joined_vertical = NULL;
        numelyx_matrix *padded = NULL;
        numelyx_matrix *zero_extended = NULL;
        numelyx_matrix *repeated = NULL;
        numelyx_matrix *block_diagonal = NULL;
        numelyx_matrix *grid_blocks[] = {NULL, NULL};
        numelyx_matrix *grid = NULL;
        numelyx_matrix *diagonal = NULL;
        numelyx_matrix *diagonal_matrix = NULL;
        numelyx_matrix *permuted = NULL;
        numelyx_matrix *vector = NULL;

        if (numelyx_matrix_create_dense(2, 2, source_values, &parent)
                != NUMELYX_SUCCESS
            || numelyx_matrix_rotate_90_clockwise(parent, &rotated)
                != NUMELYX_SUCCESS
            || numelyx_matrix_rotate_180(parent, &rotated_180)
                != NUMELYX_SUCCESS
            || numelyx_matrix_rotate_90_counterclockwise(parent, &rotated_ccw)
                != NUMELYX_SUCCESS
            || numelyx_matrix_remove_row(parent, 0, &removed)
                != NUMELYX_SUCCESS
            || numelyx_matrix_swap_columns(parent, 0, 1, &swapped)
                != NUMELYX_SUCCESS
            || numelyx_matrix_swap_rows(parent, 0, 1, &swapped_rows)
                != NUMELYX_SUCCESS
            || numelyx_matrix_permute_columns(parent, permutation, 2, &permuted_columns)
                != NUMELYX_SUCCESS
            || numelyx_matrix_join_horizontal(parent, parent, &joined)
                != NUMELYX_SUCCESS
            || numelyx_matrix_join_vertical(parent, parent, &joined_vertical)
                != NUMELYX_SUCCESS
            || numelyx_matrix_pad(parent, 1, 0, 1, 0, -1.0, &padded)
                != NUMELYX_SUCCESS
            || numelyx_matrix_zero_extend(parent, 1, 0, 1, 0, &zero_extended)
                != NUMELYX_SUCCESS
            || numelyx_matrix_repeat(parent, 2, 2, &repeated)
                != NUMELYX_SUCCESS
            || numelyx_matrix_block_diagonal(parent, 2, &block_diagonal)
                != NUMELYX_SUCCESS
            || numelyx_matrix_diagonal_extract(parent, 0, &diagonal)
                != NUMELYX_SUCCESS
            || numelyx_matrix_permute_rows(parent, permutation, 2, &permuted)
                != NUMELYX_SUCCESS
            || numelyx_matrix_create_dense(1, 2, vector_values, &vector)
                != NUMELYX_SUCCESS
            || numelyx_matrix_diagonal_from_vector(vector, 1, &diagonal_matrix)
                != NUMELYX_SUCCESS) {
            numelyx_matrix_release(diagonal_matrix);
            numelyx_matrix_release(vector);
            numelyx_matrix_release(permuted);
            numelyx_matrix_release(diagonal);
            numelyx_matrix_release(block_diagonal);
            numelyx_matrix_release(repeated);
            numelyx_matrix_release(zero_extended);
            numelyx_matrix_release(padded);
            numelyx_matrix_release(joined_vertical);
            numelyx_matrix_release(joined);
            numelyx_matrix_release(permuted_columns);
            numelyx_matrix_release(swapped_rows);
            numelyx_matrix_release(swapped);
            numelyx_matrix_release(removed);
            numelyx_matrix_release(rotated_ccw);
            numelyx_matrix_release(rotated_180);
            numelyx_matrix_release(rotated);
            numelyx_matrix_release(parent);
            return 28;
        }

        grid_blocks[0] = parent;
        grid_blocks[1] = parent;
        if (numelyx_matrix_block_grid(grid_blocks, 1, 2, &grid)
                != NUMELYX_SUCCESS) {
            numelyx_matrix_release(diagonal_matrix);
            numelyx_matrix_release(vector);
            numelyx_matrix_release(permuted);
            numelyx_matrix_release(diagonal);
            numelyx_matrix_release(block_diagonal);
            numelyx_matrix_release(repeated);
            numelyx_matrix_release(zero_extended);
            numelyx_matrix_release(padded);
            numelyx_matrix_release(joined_vertical);
            numelyx_matrix_release(joined);
            numelyx_matrix_release(permuted_columns);
            numelyx_matrix_release(swapped_rows);
            numelyx_matrix_release(swapped);
            numelyx_matrix_release(removed);
            numelyx_matrix_release(rotated_ccw);
            numelyx_matrix_release(rotated_180);
            numelyx_matrix_release(rotated);
            numelyx_matrix_release(parent);
            return 29;
        }

        numelyx_matrix_release(parent);
        numelyx_matrix_release(vector);
        if (numelyx_matrix_get(rotated, 0, 0, &value) != NUMELYX_SUCCESS
                || value != 3.0
            || numelyx_matrix_get(rotated_180, 0, 0, &value) != NUMELYX_SUCCESS
                || value != 4.0
            || numelyx_matrix_get(rotated_ccw, 0, 0, &value) != NUMELYX_SUCCESS
                || value != 2.0
            || numelyx_matrix_get(removed, 0, 0, &value) != NUMELYX_SUCCESS
                || value != 3.0
            || numelyx_matrix_get(swapped, 0, 0, &value) != NUMELYX_SUCCESS
                || value != 2.0
            || numelyx_matrix_get(swapped_rows, 0, 0, &value) != NUMELYX_SUCCESS
                || value != 3.0
            || numelyx_matrix_get(permuted_columns, 0, 0, &value)
                != NUMELYX_SUCCESS
            || value != 2.0
            || numelyx_matrix_get(joined, 0, 2, &value) != NUMELYX_SUCCESS
                || value != 1.0
            || numelyx_matrix_get(joined_vertical, 2, 1, &value)
                != NUMELYX_SUCCESS
            || value != 2.0
            || numelyx_matrix_get(padded, 0, 0, &value) != NUMELYX_SUCCESS
                || value != -1.0
            || numelyx_matrix_get(zero_extended, 0, 0, &value)
                != NUMELYX_SUCCESS
            || value != 0.0
            || numelyx_matrix_get(repeated, 2, 2, &value) != NUMELYX_SUCCESS
                || value != 1.0
            || numelyx_matrix_get(block_diagonal, 2, 2, &value)
                != NUMELYX_SUCCESS
            || value != 1.0
            || numelyx_matrix_get(grid, 1, 3, &value) != NUMELYX_SUCCESS
            || value != 4.0
            || numelyx_matrix_get(diagonal, 0, 1, &value) != NUMELYX_SUCCESS
            || value != 4.0
            || numelyx_matrix_get(diagonal_matrix, 0, 1, &value)
                != NUMELYX_SUCCESS
            || value != 5.0
            || numelyx_matrix_get(permuted, 0, 0, &value) != NUMELYX_SUCCESS
            || value != 3.0) {
            numelyx_matrix_release(diagonal_matrix);
            numelyx_matrix_release(permuted);
            numelyx_matrix_release(diagonal);
            numelyx_matrix_release(grid);
            numelyx_matrix_release(block_diagonal);
            numelyx_matrix_release(repeated);
            numelyx_matrix_release(zero_extended);
            numelyx_matrix_release(padded);
            numelyx_matrix_release(joined_vertical);
            numelyx_matrix_release(joined);
            numelyx_matrix_release(permuted_columns);
            numelyx_matrix_release(swapped_rows);
            numelyx_matrix_release(swapped);
            numelyx_matrix_release(removed);
            numelyx_matrix_release(rotated_ccw);
            numelyx_matrix_release(rotated_180);
            numelyx_matrix_release(rotated);
            return 30;
        }

        numelyx_matrix_release(diagonal_matrix);
        numelyx_matrix_release(permuted);
        numelyx_matrix_release(diagonal);
        numelyx_matrix_release(grid);
        numelyx_matrix_release(block_diagonal);
        numelyx_matrix_release(repeated);
        numelyx_matrix_release(padded);
        numelyx_matrix_release(joined);
        numelyx_matrix_release(swapped);
        numelyx_matrix_release(removed);
        numelyx_matrix_release(rotated_180);
        numelyx_matrix_release(rotated);
    }

    {
        const double coefficients[] = {4.0, 7.0, 2.0, 6.0};
        const double rhs_values[] = {1.0, 0.0};
        numelyx_matrix *coefficient_matrix = NULL;
        numelyx_matrix *rhs = NULL;
        numelyx_matrix *inverse = NULL;
        numelyx_matrix *solution = NULL;
        double determinant = 0.0;

        if (numelyx_matrix_create_dense(2, 2, coefficients, &coefficient_matrix)
                != NUMELYX_SUCCESS
            || numelyx_matrix_create_dense(2, 1, rhs_values, &rhs)
                != NUMELYX_SUCCESS
            || numelyx_matrix_determinant(coefficient_matrix, &determinant)
                != NUMELYX_SUCCESS
            || fabs(determinant - 10.0) > 1e-12
            || numelyx_matrix_inverse(coefficient_matrix, &inverse)
                != NUMELYX_SUCCESS
            || numelyx_matrix_solve(coefficient_matrix, rhs, &solution)
                != NUMELYX_SUCCESS) {
            numelyx_matrix_release(solution);
            numelyx_matrix_release(inverse);
            numelyx_matrix_release(rhs);
            numelyx_matrix_release(coefficient_matrix);
            return 17;
        }

        if (numelyx_matrix_get(inverse, 0, 0, &value) != NUMELYX_SUCCESS
                || fabs(value - 0.6) > 1e-12
            || numelyx_matrix_get(solution, 1, 0, &value) != NUMELYX_SUCCESS
                || fabs(value + 0.2) > 1e-12) {
            numelyx_matrix_release(solution);
            numelyx_matrix_release(inverse);
            numelyx_matrix_release(rhs);
            numelyx_matrix_release(coefficient_matrix);
            return 18;
        }

        numelyx_matrix_release(solution);
        numelyx_matrix_release(inverse);
        numelyx_matrix_release(rhs);
        numelyx_matrix_release(coefficient_matrix);
    }

    {
        numelyx_matrix *structured = NULL;

        if (numelyx_matrix_create_identity(3, &structured) != NUMELYX_SUCCESS
            || numelyx_matrix_get(structured, 2, 2, &value) != NUMELYX_SUCCESS
            || value != 1.0
            || numelyx_matrix_get(structured, 0, 2, &value) != NUMELYX_SUCCESS
            || value != 0.0) {
            numelyx_matrix_release(structured);
            return 9;
        }
        numelyx_matrix_release(structured);

        if (numelyx_matrix_create_zero(2, 3, &structured) != NUMELYX_SUCCESS
            || numelyx_matrix_get(structured, 1, 2, &value) != NUMELYX_SUCCESS
            || value != 0.0) {
            numelyx_matrix_release(structured);
            return 10;
        }
        numelyx_matrix_release(structured);

        if (numelyx_matrix_create_constant(2, 2, -2.5, &structured)
                != NUMELYX_SUCCESS
            || numelyx_matrix_get(structured, 1, 0, &value) != NUMELYX_SUCCESS
            || value != -2.5) {
            numelyx_matrix_release(structured);
            return 11;
        }
        numelyx_matrix_release(structured);

        if (numelyx_matrix_create_diagonal(2, (const double[]){4.0, 5.0}, &structured)
                != NUMELYX_SUCCESS
            || numelyx_matrix_get(structured, 1, 1, &value) != NUMELYX_SUCCESS
            || value != 5.0
            || numelyx_matrix_get(structured, 1, 0, &value) != NUMELYX_SUCCESS
            || value != 0.0) {
            numelyx_matrix_release(structured);
            return 12;
        }
        numelyx_matrix_release(structured);

        if (numelyx_matrix_create_scaled_identity(2, 3.0, &structured)
                != NUMELYX_SUCCESS
            || numelyx_matrix_get(structured, 1, 1, &value) != NUMELYX_SUCCESS
            || value != 3.0
            || numelyx_matrix_get(structured, 1, 0, &value) != NUMELYX_SUCCESS
            || value != 0.0) {
            numelyx_matrix_release(structured);
            return 21;
        }
        numelyx_matrix_release(structured);

        if (numelyx_matrix_create_upper_triangular(
                2,
                (const double[]){1.0, 2.0, 3.0},
                &structured
            ) != NUMELYX_SUCCESS
            || numelyx_matrix_get(structured, 0, 1, &value) != NUMELYX_SUCCESS
            || value != 2.0
            || numelyx_matrix_get(structured, 1, 0, &value) != NUMELYX_SUCCESS
            || value != 0.0) {
            numelyx_matrix_release(structured);
            return 22;
        }
        numelyx_matrix_release(structured);

        if (numelyx_matrix_create_lower_triangular(
                2,
                (const double[]){1.0, 2.0, 3.0},
                &structured
            ) != NUMELYX_SUCCESS
            || numelyx_matrix_get(structured, 1, 0, &value) != NUMELYX_SUCCESS
            || value != 2.0
            || numelyx_matrix_get(structured, 0, 1, &value) != NUMELYX_SUCCESS
            || value != 0.0) {
            numelyx_matrix_release(structured);
            return 23;
        }
        numelyx_matrix_release(structured);

        if (numelyx_matrix_create_symmetric(
                2,
                (const double[]){1.0, 2.0, 3.0},
                &structured
            ) != NUMELYX_SUCCESS
            || numelyx_matrix_get(structured, 0, 1, &value) != NUMELYX_SUCCESS
            || value != 2.0
            || numelyx_matrix_get(structured, 1, 0, &value) != NUMELYX_SUCCESS
            || value != 2.0) {
            numelyx_matrix_release(structured);
            return 24;
        }
        numelyx_matrix_release(structured);

        {
            const numelyx_sparse_entry entries[] = {
                {0, 1.0},
                {4, 5.0},
                {8, 9.0}
            };
            const numelyx_sparse_entry duplicate_entries[] = {
                {0, 1.0},
                {0, 2.0}
            };

            const double band_values[] = {
                1.0, 2.0, 0.0,
                3.0, 4.0, 5.0,
                0.0, 6.0, 7.0
            };

            if (numelyx_matrix_create_sparse(
                    3,
                    3,
                    0.0,
                    entries,
                    3,
                    &structured
                ) != NUMELYX_SUCCESS
                || numelyx_matrix_get(structured, 1, 1, &value)
                    != NUMELYX_SUCCESS
                || value != 5.0
                || numelyx_matrix_get(structured, 0, 2, &value)
                    != NUMELYX_SUCCESS
                || value != 0.0) {
                numelyx_matrix_release(structured);
                return 25;
            }
            numelyx_matrix_release(structured);
            structured = NULL;

            if (numelyx_matrix_create_banded(
                    3,
                    3,
                    1,
                    1,
                    band_values,
                    &structured
                ) != NUMELYX_SUCCESS
                || numelyx_matrix_get(structured, 1, 2, &value)
                    != NUMELYX_SUCCESS
                || value != 5.0
                || numelyx_matrix_get(structured, 0, 2, &value)
                    != NUMELYX_SUCCESS
                || value != 0.0) {
                numelyx_matrix_release(structured);
                return 27;
            }
            numelyx_matrix_release(structured);
            structured = NULL;

            if (numelyx_matrix_create_sparse(
                    2,
                    2,
                    0.0,
                    duplicate_entries,
                    2,
                    &structured
                ) != NUMELYX_INVALID_ARGUMENT
                || structured != NULL) {
                numelyx_matrix_release(structured);
                return 26;
            }
        }
    }

    numelyx_matrix_release(NULL);
    return 0;
}
