#include <numelyx/rng.h>

#include "numelyx_core_rng.h"

static numelyx_status map_rng_status(numelyx_core_rng_status status)
{
    switch (status) {
        case NUMERUS_RNG_SUCCESS:
            return NUMELYX_SUCCESS;
        case NUMERUS_RNG_INVALID_ARGUMENT:
            return NUMELYX_INVALID_ARGUMENT;
        case NUMERUS_RNG_OUT_OF_MEMORY:
            return NUMELYX_OUT_OF_MEMORY;
        case NUMERUS_RNG_SIZE_OVERFLOW:
            return NUMELYX_SIZE_OVERFLOW;
        case NUMERUS_RNG_NUMERICAL_FAILURE:
        default:
            return NUMELYX_NUMERICAL_ERROR;
    }
}

numelyx_status numelyx_rng_create(
    uint64_t seed,
    uint64_t stream,
    numelyx_rng **rng
)
{
    if (rng == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *rng = NULL;
    return map_rng_status(numelyx_core_rng_create(seed, stream, rng));
}

numelyx_status numelyx_rng_clone(
    const numelyx_rng *source,
    numelyx_rng **clone
)
{
    if (clone == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *clone = NULL;
    return map_rng_status(numelyx_core_rng_clone(source, clone));
}

numelyx_status numelyx_rng_retain(numelyx_rng *rng)
{
    return map_rng_status(numelyx_core_rng_retain(rng));
}

void numelyx_rng_release(numelyx_rng *rng)
{
    numelyx_core_rng_destroy(rng);
}

numelyx_status numelyx_rng_next_u32(numelyx_rng *rng, uint32_t *value)
{
    return map_rng_status(numelyx_core_rng_next_u32(rng, value));
}

numelyx_status numelyx_rng_uniform(numelyx_rng *rng, double *value)
{
    return map_rng_status(numelyx_core_rng_uniform(rng, value));
}

numelyx_status numelyx_rng_normal(numelyx_rng *rng, double *value)
{
    return map_rng_status(numelyx_core_rng_normal(rng, value));
}

numelyx_status numelyx_rng_sample_indices(
    numelyx_rng *rng,
    size_t population_size,
    size_t sample_size,
    bool with_replacement,
    size_t **indices
)
{
    if (indices == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *indices = NULL;
    return map_rng_status(numelyx_core_rng_sample_indices(
        rng,
        population_size,
        sample_size,
        with_replacement,
        indices
    ));
}

void numelyx_rng_free_indices(size_t *indices)
{
    numelyx_core_rng_free_indices(indices);
}
