/**
 * @file numelyx_core_distributions.h
 * @brief Internal scalar probability-distribution primitives.
 */
#ifndef NUMELYX_CORE_DISTRIBUTIONS_H
#define NUMELYX_CORE_DISTRIBUTIONS_H

#include "numelyx_core_numeric.h"

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

#endif
