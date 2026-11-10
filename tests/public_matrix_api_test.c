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

    numelyx_matrix_release(NULL);
    return 0;
}
