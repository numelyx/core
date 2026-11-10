#include <math.h>

#include <numelyx/probability.h>

int main(void)
{
    const double zero[] = {0.0};
    const double variance[] = {1.0};
    double log_density = 0.0;
    numelyx_matrix *observation = NULL;
    numelyx_matrix *mean = NULL;
    numelyx_matrix *covariance = NULL;
    numelyx_matrix *sample = NULL;
    numelyx_rng *rng = NULL;

    if (numelyx_matrix_create_dense(1, 1, zero, &observation)
            != NUMELYX_SUCCESS
        || numelyx_matrix_create_dense(1, 1, zero, &mean)
            != NUMELYX_SUCCESS
        || numelyx_matrix_create_diagonal(1, variance, &covariance)
            != NUMELYX_SUCCESS
        || numelyx_rng_create(UINT64_C(987), UINT64_C(11), &rng)
            != NUMELYX_SUCCESS) {
        numelyx_rng_release(rng);
        numelyx_matrix_release(covariance);
        numelyx_matrix_release(mean);
        numelyx_matrix_release(observation);
        return 1;
    }

    if (numelyx_multivariate_gaussian_log_density(
            observation,
            mean,
            covariance,
            &log_density
        ) != NUMELYX_SUCCESS
        || fabs(log_density + 0.9189385332046727) > 1e-12) {
        numelyx_rng_release(rng);
        numelyx_matrix_release(covariance);
        numelyx_matrix_release(mean);
        numelyx_matrix_release(observation);
        return 2;
    }

    if (numelyx_multivariate_gaussian_sample(
            mean,
            covariance,
            rng,
            &sample
        ) != NUMELYX_SUCCESS
        || sample == NULL
        || numelyx_matrix_rows(sample) != 1
        || numelyx_matrix_columns(sample) != 1) {
        numelyx_matrix_release(sample);
        numelyx_rng_release(rng);
        numelyx_matrix_release(covariance);
        numelyx_matrix_release(mean);
        numelyx_matrix_release(observation);
        return 3;
    }

    numelyx_matrix_release(sample);
    numelyx_rng_release(rng);
    numelyx_matrix_release(covariance);
    numelyx_matrix_release(mean);
    numelyx_matrix_release(observation);
    return 0;
}
