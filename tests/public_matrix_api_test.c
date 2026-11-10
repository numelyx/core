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
    }

    numelyx_matrix_release(NULL);
    return 0;
}
