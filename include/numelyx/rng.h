#ifndef NUMELYX_RNG_H
#define NUMELYX_RNG_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <numelyx/status.h>

NUMELYX_BEGIN_DECLS

/** Opaque mutable RNG state with per-object synchronized operations. */
typedef struct numelyx_core_rng numelyx_rng;

/** Create a deterministic PCG32 state from a seed and stream selector. */
NUMELYX_API numelyx_status numelyx_rng_create(
    uint64_t seed,
    uint64_t stream,
    numelyx_rng **rng
);

/** Clone the exact state and stream of source. */
NUMELYX_API numelyx_status numelyx_rng_clone(
    const numelyx_rng *source,
    numelyx_rng **clone
);

/** Retain one additional owned reference. */
NUMELYX_API numelyx_status numelyx_rng_retain(numelyx_rng *rng);

/** Release one owned reference; NULL is accepted. */
NUMELYX_API void numelyx_rng_release(numelyx_rng *rng);

/** Generate the next raw 32-bit value. */
NUMELYX_API numelyx_status numelyx_rng_next_u32(
    numelyx_rng *rng,
    uint32_t *value
);

/** Generate a uniform double in [0, 1). */
NUMELYX_API numelyx_status numelyx_rng_uniform(
    numelyx_rng *rng,
    double *value
);

/** Generate one standard-normal variate. */
NUMELYX_API numelyx_status numelyx_rng_normal(
    numelyx_rng *rng,
    double *value
);

/**
 * Sample indices from [0, population_size). The returned array belongs to the
 * caller and must be released with numelyx_rng_free_indices().
 */
NUMELYX_API numelyx_status numelyx_rng_sample_indices(
    numelyx_rng *rng,
    size_t population_size,
    size_t sample_size,
    bool with_replacement,
    size_t **indices
);

/** Release an index array returned by numelyx_rng_sample_indices(). */
NUMELYX_API void numelyx_rng_free_indices(size_t *indices);

NUMELYX_END_DECLS

#endif /* NUMELYX_RNG_H */
