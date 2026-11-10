#ifndef NUMELYX_OPTIMIZER_H
#define NUMELYX_OPTIMIZER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <numelyx/status.h>

NUMELYX_BEGIN_DECLS

typedef bool (*numelyx_objective_fn)(
    const double *parameters,
    size_t parameter_count,
    const void *context,
    double *value
);

typedef bool (*numelyx_gradient_fn)(
    const double *parameters,
    size_t parameter_count,
    const void *context,
    double *gradient
);

/* Version tag for the initial optimizer options layout. */
#define NUMELYX_OPTIMIZER_OPTIONS_VERSION_1 UINT32_C(1)

typedef struct {
    uint32_t struct_size;
    uint32_t version;
    size_t max_iterations;
    size_t max_line_search_iterations;
    double gradient_tolerance;
    double minimum_step;
    double step_tolerance;
    double initial_step;
    double finite_difference_relative_step;
} numelyx_optimizer_options;

typedef enum {
    NUMELYX_OPTIMIZER_CONVERGED_GRADIENT = 0,
    NUMELYX_OPTIMIZER_ITERATION_LIMIT,
    NUMELYX_OPTIMIZER_STEP_TOO_SMALL,
    NUMELYX_OPTIMIZER_LINE_SEARCH_FAILED,
    NUMELYX_OPTIMIZER_OBJECTIVE_CALLBACK_FAILED,
    NUMELYX_OPTIMIZER_GRADIENT_CALLBACK_FAILED,
    NUMELYX_OPTIMIZER_GRADIENT_EVALUATION_FAILED,
    NUMELYX_OPTIMIZER_NON_FINITE_OBJECTIVE,
    NUMELYX_OPTIMIZER_NON_FINITE_GRADIENT,
    NUMELYX_OPTIMIZER_NUMERICAL_BREAKDOWN
} numelyx_optimizer_termination;

typedef struct numelyx_core_optimizer_result numelyx_optimizer_result;

/** Fill default options for the BFGS minimizer. */
NUMELYX_API numelyx_status numelyx_optimizer_default_options(
    numelyx_optimizer_options *options
);

/**
 * Minimize an unconstrained smooth objective using BFGS with Armijo
 * backtracking. Initialize options with numelyx_optimizer_default_options()
 * before modifying or passing them. A NULL gradient selects central finite
 * differences. The callback context is borrowed for the duration of this call only.
 * A result is returned for valid requests even when convergence is not reached.
 */
NUMELYX_API numelyx_status numelyx_optimizer_minimize(
    numelyx_objective_fn objective,
    numelyx_gradient_fn gradient,
    const void *context,
    const double *initial_parameters,
    size_t parameter_count,
    const numelyx_optimizer_options *options,
    numelyx_optimizer_result **result
);

NUMELYX_API numelyx_status numelyx_optimizer_result_retain(
    numelyx_optimizer_result *result
);
NUMELYX_API void numelyx_optimizer_result_release(
    numelyx_optimizer_result *result
);
NUMELYX_API size_t numelyx_optimizer_result_parameter_count(
    const numelyx_optimizer_result *result
);
NUMELYX_API numelyx_status numelyx_optimizer_result_get_parameter(
    const numelyx_optimizer_result *result,
    size_t index,
    double *value
);
NUMELYX_API bool numelyx_optimizer_result_objective_value_valid(
    const numelyx_optimizer_result *result
);
NUMELYX_API numelyx_status numelyx_optimizer_result_get_objective(
    const numelyx_optimizer_result *result,
    double *value
);
NUMELYX_API numelyx_status numelyx_optimizer_result_get_termination(
    const numelyx_optimizer_result *result,
    numelyx_optimizer_termination *termination
);
NUMELYX_API size_t numelyx_optimizer_result_iterations(
    const numelyx_optimizer_result *result
);
NUMELYX_API size_t numelyx_optimizer_result_objective_evaluations(
    const numelyx_optimizer_result *result
);
NUMELYX_API size_t numelyx_optimizer_result_gradient_evaluations(
    const numelyx_optimizer_result *result
);

NUMELYX_END_DECLS

#endif /* NUMELYX_OPTIMIZER_H */
