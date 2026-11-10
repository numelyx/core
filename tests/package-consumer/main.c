#include <numelyx/numelyx.h>

int main(void)
{
    const double values[] = {2.0, 4.0, 6.0, 8.0};
    double element = 0.0;
    numelyx_matrix *matrix = NULL;
    numelyx_rng *rng = NULL;

    if (numelyx_core_version_major() != (uint32_t) NUMELYX_EXPECTED_VERSION_MAJOR
        || numelyx_core_version_minor() != (uint32_t) NUMELYX_EXPECTED_VERSION_MINOR
        || numelyx_core_version_patch() != (uint32_t) NUMELYX_EXPECTED_VERSION_PATCH
        || numelyx_core_abi_version() != UINT32_C(0)) {
        return 1;
    }

    if (numelyx_matrix_create_dense(2, 2, values, &matrix) != NUMELYX_SUCCESS) {
        return 2;
    }
    if (numelyx_matrix_get(matrix, 1, 0, &element) != NUMELYX_SUCCESS
        || element != 6.0) {
        numelyx_matrix_release(matrix);
        return 3;
    }
    numelyx_matrix_release(matrix);

    if (numelyx_rng_create(UINT64_C(42), UINT64_C(54), &rng)
        != NUMELYX_SUCCESS) {
        return 4;
    }

    numelyx_rng_release(rng);
    return 0;
}
