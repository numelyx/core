/**
 * @file numelyx_core_distributions.h
 * @brief Internal scalar probability-distribution primitives.
 */
#ifndef NUMELYX_CORE_DISTRIBUTIONS_H
#define NUMELYX_CORE_DISTRIBUTIONS_H

#include "numelyx_core_numeric.h"
#include "numelyx_core_rng.h"

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

#endif
