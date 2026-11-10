/**
 * @file numelyx_core_distributions.h
 * @brief Internal scalar probability-distribution primitives.
 */
#ifndef NUMELYX_CORE_DISTRIBUTIONS_H
#define NUMELYX_CORE_DISTRIBUTIONS_H

#include "numelyx_core_numeric.h"
#include "numelyx_core_rng.h"
#include <stdint.h>

/** Normal distribution with finite location and positive finite scale. */
numelyx_core_numeric_status numelyx_core_normal_log_pdf(
    double x, double location, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_normal_pdf(
    double x, double location, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_normal_log_cdf(
    double x, double location, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_normal_log_survival(
    double x, double location, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_normal_cdf(
    double x, double location, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_normal_survival(
    double x, double location, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_normal_distribution_quantile(
    double probability, double location, double scale, double *result
);

/** Student-t distribution with finite location, positive scale and degrees of freedom. */
numelyx_core_numeric_status numelyx_core_student_t_log_pdf(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
);
numelyx_core_numeric_status numelyx_core_student_t_pdf(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
);
numelyx_core_numeric_status numelyx_core_student_t_log_cdf(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
);
numelyx_core_numeric_status numelyx_core_student_t_log_survival(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
);
numelyx_core_numeric_status numelyx_core_student_t_cdf(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
);
numelyx_core_numeric_status numelyx_core_student_t_survival(
    double x, double degrees_of_freedom, double location, double scale,
    double *result
);
numelyx_core_numeric_status numelyx_core_student_t_quantile(
    double probability, double degrees_of_freedom, double location,
    double scale, double *result
);

/** Chi-square distribution with finite positive degrees of freedom. */
numelyx_core_numeric_status numelyx_core_chi_square_log_pdf(
    double x, double degrees_of_freedom, double *result
);
numelyx_core_numeric_status numelyx_core_chi_square_pdf(
    double x, double degrees_of_freedom, double *result
);
numelyx_core_numeric_status numelyx_core_chi_square_log_cdf(
    double x, double degrees_of_freedom, double *result
);
numelyx_core_numeric_status numelyx_core_chi_square_log_survival(
    double x, double degrees_of_freedom, double *result
);
numelyx_core_numeric_status numelyx_core_chi_square_cdf(
    double x, double degrees_of_freedom, double *result
);
numelyx_core_numeric_status numelyx_core_chi_square_survival(
    double x, double degrees_of_freedom, double *result
);
numelyx_core_numeric_status numelyx_core_chi_square_quantile(
    double probability, double degrees_of_freedom, double *result
);

/** Draw a normal variate from an explicit RNG; output and RNG are unchanged on failure. */
numelyx_core_numeric_status numelyx_core_normal_sample(
    numelyx_core_rng *rng, double location, double scale, double *result
);

/** Gamma distribution parameterized by positive shape and scale. */
numelyx_core_numeric_status numelyx_core_gamma_log_pdf(
    double x, double shape, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_gamma_pdf(
    double x, double shape, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_gamma_log_cdf(
    double x, double shape, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_gamma_log_survival(
    double x, double shape, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_gamma_cdf(
    double x, double shape, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_gamma_survival(
    double x, double shape, double scale, double *result
);
numelyx_core_numeric_status numelyx_core_gamma_quantile(
    double probability, double shape, double scale, double *result
);

/** Beta distribution on [0,1] with positive alpha and beta parameters. */
numelyx_core_numeric_status numelyx_core_beta_log_pdf(
    double x, double alpha, double beta, double *result
);
numelyx_core_numeric_status numelyx_core_beta_pdf(
    double x, double alpha, double beta, double *result
);
numelyx_core_numeric_status numelyx_core_beta_log_cdf(
    double x, double alpha, double beta, double *result
);
numelyx_core_numeric_status numelyx_core_beta_log_survival(
    double x, double alpha, double beta, double *result
);
numelyx_core_numeric_status numelyx_core_beta_cdf(
    double x, double alpha, double beta, double *result
);
numelyx_core_numeric_status numelyx_core_beta_survival(
    double x, double alpha, double beta, double *result
);
numelyx_core_numeric_status numelyx_core_beta_quantile(
    double probability, double alpha, double beta, double *result
);

/** F distribution with positive numerator and denominator degrees of freedom. */
numelyx_core_numeric_status numelyx_core_f_log_pdf(
    double x, double numerator_df, double denominator_df, double *result
);
numelyx_core_numeric_status numelyx_core_f_pdf(
    double x, double numerator_df, double denominator_df, double *result
);
numelyx_core_numeric_status numelyx_core_f_log_cdf(
    double x, double numerator_df, double denominator_df, double *result
);
numelyx_core_numeric_status numelyx_core_f_log_survival(
    double x, double numerator_df, double denominator_df, double *result
);
numelyx_core_numeric_status numelyx_core_f_cdf(
    double x, double numerator_df, double denominator_df, double *result
);
numelyx_core_numeric_status numelyx_core_f_survival(
    double x, double numerator_df, double denominator_df, double *result
);
numelyx_core_numeric_status numelyx_core_f_quantile(
    double probability, double numerator_df, double denominator_df,
    double *result
);

/** Discrete Poisson distribution; counts above 2^53 are not supported. */
numelyx_core_numeric_status numelyx_core_poisson_log_pmf(
    uint64_t count, double rate, double *result
);
numelyx_core_numeric_status numelyx_core_poisson_pmf(
    uint64_t count, double rate, double *result
);
numelyx_core_numeric_status numelyx_core_poisson_log_cdf(
    uint64_t count, double rate, double *result
);
numelyx_core_numeric_status numelyx_core_poisson_log_survival(
    uint64_t count, double rate, double *result
);
numelyx_core_numeric_status numelyx_core_poisson_cdf(
    uint64_t count, double rate, double *result
);
numelyx_core_numeric_status numelyx_core_poisson_survival(
    uint64_t count, double rate, double *result
);
numelyx_core_numeric_status numelyx_core_poisson_quantile(
    double probability, double rate, uint64_t *result
);

/** Binomial distribution; trial counts above 2^53 are not supported. */
numelyx_core_numeric_status numelyx_core_binomial_log_pmf(
    uint64_t successes, uint64_t trials, double probability, double *result
);
numelyx_core_numeric_status numelyx_core_binomial_pmf(
    uint64_t successes, uint64_t trials, double probability, double *result
);
numelyx_core_numeric_status numelyx_core_binomial_log_cdf(
    uint64_t successes, uint64_t trials, double probability, double *result
);
numelyx_core_numeric_status numelyx_core_binomial_log_survival(
    uint64_t successes, uint64_t trials, double probability, double *result
);
numelyx_core_numeric_status numelyx_core_binomial_cdf(
    uint64_t successes, uint64_t trials, double probability, double *result
);
numelyx_core_numeric_status numelyx_core_binomial_survival(
    uint64_t successes, uint64_t trials, double probability, double *result
);
numelyx_core_numeric_status numelyx_core_binomial_quantile(
    double quantile_probability, uint64_t trials, double success_probability,
    uint64_t *result
);

#endif
