/**
 * @file numelyx_core_root.h
 * @brief Allocation-free scalar root-finding primitives.
 */
#ifndef NUMELYX_CORE_ROOT_H
#define NUMELYX_CORE_ROOT_H

#include <stdbool.h>
#include <stddef.h>

/** Callback returns true only when it writes a complete function value. */
typedef bool (*numelyx_core_scalar_function_fn)(
    double x, const void *context, double *value
);

typedef enum {
    NUMELYX_ROOT_SUCCESS = 0,
    NUMELYX_ROOT_INVALID_ARGUMENT,
    NUMELYX_ROOT_NON_FINITE_INPUT,
    NUMELYX_ROOT_NO_BRACKET,
    NUMELYX_ROOT_CALLBACK_FAILED,
    NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE,
    NUMELYX_ROOT_NON_FINITE_DERIVATIVE,
    NUMELYX_ROOT_ITERATION_LIMIT,
    NUMELYX_ROOT_NUMERICAL_FAILURE
} numelyx_core_root_status;

typedef bool (*numelyx_core_scalar_derivative_fn)(
    double x, const void *context, double *derivative
);

typedef struct {
    size_t max_iterations;
    double absolute_tolerance;
    double relative_tolerance;
} numelyx_core_bisection_options;

typedef struct {
    double root;
    double function_value;
    size_t iterations;
    size_t function_evaluations;
    size_t derivative_evaluations;
} numelyx_core_root_result;

/** Fill conservative defaults; output is unchanged when options is NULL. */
numelyx_core_root_status numelyx_core_bisection_default_options(
    numelyx_core_bisection_options *options
);

/**
 * Find a root in [lower, upper] using bracket-preserving bisection.
 * Function values at both endpoints must have opposite signs or include zero.
 * On iteration-limit termination, result contains the best midpoint found.
 * No allocation is performed.
 */
numelyx_core_root_status numelyx_core_bisect(
    numelyx_core_scalar_function_fn function,
    const void *context,
    double lower,
    double upper,
    const numelyx_core_bisection_options *options,
    numelyx_core_root_result *result
);

numelyx_core_root_status numelyx_core_safeguarded_newton(
    numelyx_core_scalar_function_fn function,
    numelyx_core_scalar_derivative_fn derivative,
    const void *context,
    double lower,
    double upper,
    double initial,
    const numelyx_core_bisection_options *options,
    numelyx_core_root_result *result
);

#endif
