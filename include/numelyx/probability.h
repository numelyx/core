#ifndef NUMELYX_PROBABILITY_H
#define NUMELYX_PROBABILITY_H

#include <numelyx/matrix.h>
#include <numelyx/rng.h>

NUMELYX_BEGIN_DECLS

/**
 * Compute multivariate Gaussian log density for column-vector observation and
 * mean with symmetric positive-definite covariance. Uses Cholesky and log
 * determinant; it does not form the covariance inverse. Inputs are borrowed.
 */
NUMELYX_API numelyx_status numelyx_multivariate_gaussian_log_density(
    const numelyx_matrix *observation,
    const numelyx_matrix *mean,
    const numelyx_matrix *covariance,
    double *log_density
);

/**
 * Draw one multivariate Gaussian sample using the supplied RNG state.
 * The returned Matrix is owned by the caller and must be released.
 */
NUMELYX_API numelyx_status numelyx_multivariate_gaussian_sample(
    const numelyx_matrix *mean,
    const numelyx_matrix *covariance,
    numelyx_rng *rng,
    numelyx_matrix **sample
);

NUMELYX_END_DECLS

#endif /* NUMELYX_PROBABILITY_H */
