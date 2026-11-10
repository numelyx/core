#include "../src/core/numelyx_core_root.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

static bool square_minus_two(double x, const void *context, double *value)
{
    (void) context;
    *value = x * x - 2.0;
    return true;
}

static bool linear_root(double x, const void *context, double *value)
{
    double root = *(const double *) context;
    *value = x - root;
    return true;
}

static bool callback_failure(double x, const void *context, double *value)
{
    (void) x;
    (void) context;
    (void) value;
    return false;
}

static bool nonfinite_value(double x, const void *context, double *value)
{
    (void) x;
    (void) context;
    *value = NAN;
    return true;
}

static bool same_sign(double x, const void *context, double *value)
{
    (void) context;
    *value = x * x + 1.0;
    return true;
}

static void test_bisection_reference_and_diagnostics(void)
{
    numelyx_core_bisection_options options;
    numelyx_core_root_result result;
    assert(numelyx_core_bisection_default_options(&options) ==
        NUMELYX_ROOT_SUCCESS);
    assert(numelyx_core_bisect(
        square_minus_two, NULL, 0.0, 2.0, &options, &result) ==
        NUMELYX_ROOT_SUCCESS);
    assert(fabs(result.root - sqrt(2.0)) < 2e-12);
    assert(fabs(result.function_value) < 6e-12);
    assert(result.iterations > 0);
    assert(result.function_evaluations == result.iterations + 2);
}

static void test_endpoint_roots_and_extreme_brackets(void)
{
    numelyx_core_bisection_options options;
    numelyx_core_root_result result;
    double root = -DBL_MAX * 0.5;
    assert(numelyx_core_bisection_default_options(&options) ==
        NUMELYX_ROOT_SUCCESS);
    assert(numelyx_core_bisect(
        linear_root, &root, root, DBL_MAX, &options, &result) ==
        NUMELYX_ROOT_SUCCESS);
    assert(result.root == root);
    assert(result.iterations == 0);
    root = 0.0;
    assert(numelyx_core_bisect(
        linear_root, &root, -DBL_MAX, DBL_MAX, &options, &result) ==
        NUMELYX_ROOT_SUCCESS);
    assert(result.root == 0.0);
}

static void test_bisection_failure_contracts(void)
{
    numelyx_core_bisection_options options;
    numelyx_core_root_result result = {11.0, 12.0, 13, 14, 15};
    assert(numelyx_core_bisection_default_options(&options) ==
        NUMELYX_ROOT_SUCCESS);

    assert(numelyx_core_bisect(
        same_sign, NULL, -1.0, 1.0, &options, &result) ==
        NUMELYX_ROOT_NO_BRACKET);
    assert(result.root == 11.0 && result.function_value == 12.0);

    assert(numelyx_core_bisect(
        callback_failure, NULL, -1.0, 1.0, &options, &result) ==
        NUMELYX_ROOT_CALLBACK_FAILED);
    assert(result.root == 11.0);

    assert(numelyx_core_bisect(
        nonfinite_value, NULL, -1.0, 1.0, &options, &result) ==
        NUMELYX_ROOT_NON_FINITE_FUNCTION_VALUE);
    assert(result.root == 11.0);

    options.max_iterations = 1;
    options.absolute_tolerance = DBL_MIN;
    options.relative_tolerance = 0.0;
    assert(numelyx_core_bisect(
        square_minus_two, NULL, 0.0, 2.0, &options, &result) ==
        NUMELYX_ROOT_ITERATION_LIMIT);
    assert(result.iterations == 1);
    assert(result.function_evaluations == 4);

    assert(numelyx_core_bisect(
        square_minus_two, NULL, 1.0, 1.0, &options, &result) ==
        NUMELYX_ROOT_INVALID_ARGUMENT);
    assert(numelyx_core_bisection_default_options(NULL) ==
        NUMELYX_ROOT_INVALID_ARGUMENT);
    assert(numelyx_core_bisect(
        square_minus_two, NULL, -INFINITY, 2.0, &options, &result) ==
        NUMELYX_ROOT_NON_FINITE_INPUT);
    assert(result.root == 11.0);
    options.absolute_tolerance = NAN;
    assert(numelyx_core_bisect(
        square_minus_two, NULL, 0.0, 2.0, &options, &result) ==
        NUMELYX_ROOT_NON_FINITE_INPUT);
    assert(result.root == 11.0);
}


static bool square_derivative(double x, const void *context, double *value)
{
    (void) context;
    *value = 2.0 * x;
    return true;
}

static bool zero_derivative(double x, const void *context, double *value)
{
    (void) x;
    (void) context;
    *value = 0.0;
    return true;
}

static bool failing_derivative(double x, const void *context, double *value)
{
    (void) x;
    (void) context;
    (void) value;
    return false;
}

static bool nonfinite_derivative(double x, const void *context, double *value)
{
    (void) x;
    (void) context;
    *value = NAN;
    return true;
}

static void test_safeguarded_newton_and_bisection_fallback(void)
{
    numelyx_core_bisection_options options;
    numelyx_core_root_result result;
    assert(numelyx_core_bisection_default_options(&options) ==
        NUMELYX_ROOT_SUCCESS);
    assert(numelyx_core_safeguarded_newton(
        square_minus_two, square_derivative, NULL,
        0.0, 2.0, 1.0, &options, &result) == NUMELYX_ROOT_SUCCESS);
    assert(fabs(result.root - sqrt(2.0)) < 2e-12);
    assert(result.derivative_evaluations > 0);
    assert(result.function_evaluations >= 3);

    assert(numelyx_core_safeguarded_newton(
        square_minus_two, zero_derivative, NULL,
        0.0, 2.0, 1.0, &options, &result) == NUMELYX_ROOT_SUCCESS);
    assert(fabs(result.root - sqrt(2.0)) < 2e-12);
    assert(result.derivative_evaluations > 0);
}

static void test_safeguarded_newton_failure_contracts(void)
{
    numelyx_core_bisection_options options;
    numelyx_core_root_result result = {1.0, 2.0, 3, 4, 5};
    assert(numelyx_core_bisection_default_options(&options) ==
        NUMELYX_ROOT_SUCCESS);
    assert(numelyx_core_safeguarded_newton(
        square_minus_two, failing_derivative, NULL,
        0.0, 2.0, 1.0, &options, &result) ==
        NUMELYX_ROOT_CALLBACK_FAILED);
    assert(result.root == 1.0 && result.derivative_evaluations == 5);
    assert(numelyx_core_safeguarded_newton(
        square_minus_two, nonfinite_derivative, NULL,
        0.0, 2.0, 1.0, &options, &result) ==
        NUMELYX_ROOT_NON_FINITE_DERIVATIVE);
    assert(result.root == 1.0);
    assert(numelyx_core_safeguarded_newton(
        square_minus_two, square_derivative, NULL,
        0.0, 2.0, 3.0, &options, &result) ==
        NUMELYX_ROOT_INVALID_ARGUMENT);
    assert(numelyx_core_safeguarded_newton(
        square_minus_two, square_derivative, NULL,
        0.0, 2.0, NAN, &options, &result) ==
        NUMELYX_ROOT_NON_FINITE_INPUT);
}

int main(void)
{
    test_bisection_reference_and_diagnostics();
    test_endpoint_roots_and_extreme_brackets();
    test_bisection_failure_contracts();
    test_safeguarded_newton_and_bisection_fallback();
    test_safeguarded_newton_failure_contracts();
    puts("Root-finding tests passed.");
    return 0;
}
