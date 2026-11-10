/**
 * @file numelyx_core_root.c
 * @brief Robust scalar bisection with explicit convergence diagnostics.
 */
#include "numelyx_core_root.h"

#include <float.h>
#include <math.h>

static bool valid_options(const numelyx_core_bisection_options *options)
{
    return options != NULL &&
        options->max_iterations > 0 &&
        isfinite(options->absolute_tolerance) &&
        isfinite(options->relative_tolerance) &&
        options->absolute_tolerance >= 0.0 &&
        options->relative_tolerance >= 0.0 &&
        (options->absolute_tolerance > 0.0 ||
         options->relative_tolerance > 0.0);
}

static double safe_midpoint(double lower, double upper)
{
    if (signbit(lower) != signbit(upper)) {
        return lower * 0.5 + upper * 0.5;
    }
    return lower + (upper - lower) * 0.5;
}

static bool interval_converged(
    double lower, double upper,
    const numelyx_core_bisection_options *options
)
{
    double half_width = fabs(upper * 0.5 - lower * 0.5);
    double scale = fmax(fabs(lower), fabs(upper));
    double threshold = options->absolute_tolerance * 0.5 +
        options->relative_tolerance * (scale * 0.5);
    return half_width <= threshold;
}

numelyx_core_root_status numelyx_core_bisection_default_options(
    numelyx_core_bisection_options *options
)
{
    if (options == NULL) return NUMELYX_ROOT_INVALID_ARGUMENT;
    options->max_iterations = 256;
    options->absolute_tolerance = 1e-12;
    options->relative_tolerance = 8.0 * DBL_EPSILON;
    return NUMELYX_ROOT_SUCCESS;
}

numelyx_core_root_status numelyx_core_bisect(
    numelyx_core_scalar_function_fn function,
    const void *context,
    double lower,
    double upper,
    const numelyx_core_bisection_options *options,
    numelyx_core_root_result *result
)
{
    double lower_value;
    double upper_value;
    size_t evaluations = 0;
    size_t iteration;

    if (function == NULL || result == NULL || !valid_options(options) ||
        !isfinite(lower) || !isfinite(upper) || !(lower < upper)) {
        return NUMELYX_ROOT_INVALID_ARGUMENT;
    }

    if (!function(lower, context, &lower_value)) {
        return NUMELYX_ROOT_CALLBACK_FAILED;
    }
    evaluations++;
    if (!isfinite(lower_value)) return NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE;
    if (lower_value == 0.0) {
        result->root = lower;
        result->function_value = lower_value;
        result->iterations = 0;
        result->function_evaluations = evaluations;
        result->derivative_evaluations = 0;
        return NUMELYX_ROOT_SUCCESS;
    }

    if (!function(upper, context, &upper_value)) {
        return NUMELYX_ROOT_CALLBACK_FAILED;
    }
    evaluations++;
    if (!isfinite(upper_value)) return NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE;
    if (upper_value == 0.0) {
        result->root = upper;
        result->function_value = upper_value;
        result->iterations = 0;
        result->function_evaluations = evaluations;
        return NUMELYX_ROOT_SUCCESS;
    }

    if (signbit(lower_value) == signbit(upper_value))
        return NUMELYX_ROOT_NO_BRACKET;

    for (iteration = 1; iteration <= options->max_iterations; iteration++) {
        double midpoint = safe_midpoint(lower, upper);
        double midpoint_value;

        if (midpoint == lower || midpoint == upper) {
            double best = fabs(lower_value) <= fabs(upper_value) ? lower : upper;
            double best_value = best == lower ? lower_value : upper_value;
            result->root = best;
            result->function_value = best_value;
            result->iterations = iteration - 1;
            result->function_evaluations = evaluations;
            result->derivative_evaluations = 0;
            return NUMELYX_ROOT_SUCCESS;
        }

        if (!function(midpoint, context, &midpoint_value)) {
            return NUMELYX_ROOT_CALLBACK_FAILED;
        }
        evaluations++;
        if (!isfinite(midpoint_value))
            return NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE;

        if (midpoint_value == 0.0 || interval_converged(lower, upper, options)) {
            result->root = midpoint;
            result->function_value = midpoint_value;
            result->iterations = iteration;
            result->function_evaluations = evaluations;
            return NUMELYX_ROOT_SUCCESS;
        }

        if (signbit(midpoint_value) == signbit(lower_value)) {
            lower = midpoint;
            lower_value = midpoint_value;
        } else {
            upper = midpoint;
            upper_value = midpoint_value;
        }
    }

    {
        double midpoint = safe_midpoint(lower, upper);
        double midpoint_value;
        if (!function(midpoint, context, &midpoint_value)) {
            return NUMELYX_ROOT_CALLBACK_FAILED;
        }
        evaluations++;
        if (!isfinite(midpoint_value))
            return NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE;
        result->root = midpoint;
        result->function_value = midpoint_value;
        result->iterations = options->max_iterations;
        result->function_evaluations = evaluations;
        result->derivative_evaluations = 0;
    }
    return NUMELYX_ROOT_ITERATION_LIMIT;
}

 
numelyx_core_root_status numelyx_core_safeguarded_newton(
    numelyx_core_scalar_function_fn function,
    numelyx_core_scalar_derivative_fn derivative,
    const void *context,
    double lower,
    double upper,
    double initial,
    const numelyx_core_bisection_options *options,
    numelyx_core_root_result *result
)
{
    double lower_value;
    double upper_value;
    double x;
    double fx;
    size_t function_evaluations = 0;
    size_t derivative_evaluations = 0;
    size_t iteration;

    if (function == NULL || derivative == NULL || result == NULL ||
        !valid_options(options) || !isfinite(lower) || !isfinite(upper) ||
        !isfinite(initial) || !(lower < upper) ||
        !(initial > lower && initial < upper)) {
        return NUMELYX_ROOT_INVALID_ARGUMENT;
    }

    if (!function(lower, context, &lower_value))
        return NUMELYX_ROOT_CALLBACK_FAILED;
    function_evaluations++;
    if (!isfinite(lower_value))
        return NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE;
    if (lower_value == 0.0) {
        result->root = lower;
        result->function_value = lower_value;
        result->iterations = 0;
        result->function_evaluations = function_evaluations;
        result->derivative_evaluations = 0;
        return NUMELYX_ROOT_SUCCESS;
    }

    if (!function(upper, context, &upper_value))
        return NUMELYX_ROOT_CALLBACK_FAILED;
    function_evaluations++;
    if (!isfinite(upper_value))
        return NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE;
    if (upper_value == 0.0) {
        result->root = upper;
        result->function_value = upper_value;
        result->iterations = 0;
        result->function_evaluations = function_evaluations;
        result->derivative_evaluations = 0;
        return NUMELYX_ROOT_SUCCESS;
    }
    if (signbit(lower_value) == signbit(upper_value))
        return NUMELYX_ROOT_NO_BRACKET;

    x = initial;
    if (!function(x, context, &fx))
        return NUMELYX_ROOT_CALLBACK_FAILED;
    function_evaluations++;
    if (!isfinite(fx))
        return NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE;
    if (fx == 0.0) {
        result->root = x;
        result->function_value = fx;
        result->iterations = 0;
        result->function_evaluations = function_evaluations;
        result->derivative_evaluations = 0;
        return NUMELYX_ROOT_SUCCESS;
    }

    for (iteration = 1; iteration <= options->max_iterations; iteration++) {
        double slope;
        double candidate;
        double candidate_value;

        if (interval_converged(lower, upper, options)) {
            result->root = x;
            result->function_value = fx;
            result->iterations = iteration - 1;
            result->function_evaluations = function_evaluations;
            result->derivative_evaluations = derivative_evaluations;
            return NUMELYX_ROOT_SUCCESS;
        }

        derivative_evaluations++;
        if (!derivative(x, context, &slope))
            return NUMELYX_ROOT_CALLBACK_FAILED;
        if (!isfinite(slope))
            return NUMELYX_ROOT_NON_FINITE_DERIVATIVE;

        candidate = slope == 0.0 ? NAN : x - fx / slope;
        if (!isfinite(candidate) || candidate <= lower || candidate >= upper ||
            candidate == x) {
            candidate = safe_midpoint(lower, upper);
        }
        if (candidate == lower || candidate == upper) {
            result->root = x;
            result->function_value = fx;
            result->iterations = iteration - 1;
            result->function_evaluations = function_evaluations;
            result->derivative_evaluations = derivative_evaluations;
            return NUMELYX_ROOT_SUCCESS;
        }

        if (!function(candidate, context, &candidate_value))
            return NUMELYX_ROOT_CALLBACK_FAILED;
        function_evaluations++;
        if (!isfinite(candidate_value))
            return NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE;

        if (candidate_value == 0.0) {
            result->root = candidate;
            result->function_value = candidate_value;
            result->iterations = iteration;
            result->function_evaluations = function_evaluations;
            result->derivative_evaluations = derivative_evaluations;
            return NUMELYX_ROOT_SUCCESS;
        }

        if (signbit(candidate_value) == signbit(lower_value)) {
            lower = candidate;
            lower_value = candidate_value;
        } else {
            upper = candidate;
            upper_value = candidate_value;
        }
        x = candidate;
        fx = candidate_value;
    }

    result->root = x;
    result->function_value = fx;
    result->iterations = options->max_iterations;
    result->function_evaluations = function_evaluations;
    result->derivative_evaluations = derivative_evaluations;
    return NUMELYX_ROOT_ITERATION_LIMIT;
}
