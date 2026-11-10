#include <numelyx/rng.h>

int main(void)
{
    numelyx_rng *first = NULL;
    numelyx_rng *clone = NULL;
    uint32_t first_value = 0;
    uint32_t clone_value = 0;
    double uniform = -1.0;
    size_t *indices = NULL;

    if (numelyx_rng_create(UINT64_C(1234), UINT64_C(7), &first)
            != NUMELYX_SUCCESS
        || numelyx_rng_clone(first, &clone) != NUMELYX_SUCCESS) {
        numelyx_rng_release(clone);
        numelyx_rng_release(first);
        return 1;
    }

    if (numelyx_rng_next_u32(first, &first_value) != NUMELYX_SUCCESS
        || numelyx_rng_next_u32(clone, &clone_value) != NUMELYX_SUCCESS
        || first_value != clone_value) {
        numelyx_rng_release(clone);
        numelyx_rng_release(first);
        return 2;
    }

    if (numelyx_rng_uniform(first, &uniform) != NUMELYX_SUCCESS
        || !(uniform >= 0.0 && uniform < 1.0)) {
        numelyx_rng_release(clone);
        numelyx_rng_release(first);
        return 3;
    }

    if (numelyx_rng_sample_indices(first, 20, 5, false, &indices)
            != NUMELYX_SUCCESS
        || indices == NULL) {
        numelyx_rng_free_indices(indices);
        numelyx_rng_release(clone);
        numelyx_rng_release(first);
        return 4;
    }
    for (size_t i = 0; i < 5; i++) {
        if (indices[i] >= 20) {
            numelyx_rng_free_indices(indices);
            numelyx_rng_release(clone);
            numelyx_rng_release(first);
            return 5;
        }
        for (size_t j = 0; j < i; j++) {
            if (indices[i] == indices[j]) {
                numelyx_rng_free_indices(indices);
                numelyx_rng_release(clone);
                numelyx_rng_release(first);
                return 6;
            }
        }
    }
    numelyx_rng_free_indices(indices);

    numelyx_rng_release(clone);
    numelyx_rng_release(first);
    numelyx_rng_release(NULL);
    return 0;
}
