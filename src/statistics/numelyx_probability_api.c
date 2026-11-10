#include <numelyx/probability.h>

#include "numelyx_core_probability.h"

static numelyx_status map_probability_status(
    numelyx_core_probability_status status
)
{
    switch (status) {
        case NUMERUS_PROBABILITY_SUCCESS:
            return NUMELYX_SUCCESS;
        case NUMERUS_PROBABILITY_INVALID_ARGUMENT:
            return NUMELYX_INVALID_ARGUMENT;
        case NUMERUS_PROBABILITY_DIMENSION_MISMATCH:
            return NUMELYX_DIMENSION_MISMATCH;
        case NUMERUS_PROBABILITY_NON_FINITE_INPUT:
            return NUMELYX_NON_FINITE;
        case NUMERUS_PROBABILITY_NOT_POSITIVE_DEFINITE:
            return NUMELYX_NOT_POSITIVE_DEFINITE;
        case NUMERUS_PROBABILITY_SIZE_OVERFLOW:
            return NUMELYX_SIZE_OVERFLOW;
        case NUMERUS_PROBABILITY_OUT_OF_MEMORY:
            return NUMELYX_OUT_OF_MEMORY;
        case NUMERUS_PROBABILITY_NUMERICAL_FAILURE:
        default:
            return NUMELYX_NUMERICAL_ERROR;
    }
}

numelyx_status numelyx_multivariate_gaussian_log_density(
    const numelyx_matrix *observation,
    const numelyx_matrix *mean,
    const numelyx_matrix *covariance,
    double *log_density
)
{
    return map_probability_status(
        numelyx_core_multivariate_gaussian_log_density(
            observation,
            mean,
            covariance,
            log_density
        )
    );
}

numelyx_status numelyx_multivariate_gaussian_sample(
    const numelyx_matrix *mean,
    const numelyx_matrix *covariance,
    numelyx_rng *rng,
    numelyx_matrix **sample
)
{
    if (sample == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *sample = NULL;

    return map_probability_status(numelyx_core_multivariate_gaussian_sample(
        mean,
        covariance,
        rng,
        sample
    ));
}
