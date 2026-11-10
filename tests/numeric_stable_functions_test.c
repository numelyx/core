/**
 * @file numeric_stable_functions_test.c
 * @brief Native tests for stable sigmoid/logit/softplus and log-domain numerical functions.
 *
 * @details These native tests define regression coverage for the named API
 * contract. Assertions are executable specifications: intentional behavior
 * changes should update these checks together with the corresponding API
 * documentation. This file is test support, not runtime code.
 */

#include "../numelyx_core_numeric.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

static void assert_close(double actual, double expected, double tolerance)
{
    assert(fabs(actual - expected) <= tolerance);
}

static void test_sigmoid_stability(void)
{
    double value = 123.0;

    assert(numelyx_core_numeric_sigmoid(0.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.5);
    assert(numelyx_core_numeric_sigmoid(DBL_MAX, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == 1.0);
    assert(numelyx_core_numeric_sigmoid(-DBL_MAX, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.0);
    assert(numelyx_core_numeric_sigmoid(NAN, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 0.0);
    assert(numelyx_core_numeric_sigmoid(INFINITY, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 0.0);
    assert(numelyx_core_numeric_sigmoid(0.0, NULL) ==
        NUMERUS_NUMERIC_INVALID_ARGUMENT);
}

static void test_logit_domain_and_endpoints(void)
{
    double value = 77.0;
    double near_zero = nextafter(0.0, 1.0);
    double near_one = nextafter(1.0, 0.0);

    assert(numelyx_core_numeric_logit(0.5, &value) == NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.0);
    assert(numelyx_core_numeric_logit(near_zero, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(isfinite(value) && value < 0.0);
    assert(numelyx_core_numeric_logit(near_one, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(isfinite(value) && value > 0.0);

    value = 77.0;
    assert(numelyx_core_numeric_logit(0.0, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 77.0);
    assert(numelyx_core_numeric_logit(1.0, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 77.0);
    assert(numelyx_core_numeric_logit(-0.1, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 77.0);
    assert(numelyx_core_numeric_logit(1.1, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 77.0);
    assert(numelyx_core_numeric_logit(NAN, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 77.0);
    assert(numelyx_core_numeric_logit(0.5, NULL) ==
        NUMERUS_NUMERIC_INVALID_ARGUMENT);
}

static void test_softplus_and_log_sigmoid_extremes(void)
{
    double value = 123.0;

    assert(numelyx_core_numeric_softplus(0.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert_close(value, log(2.0), 1e-15);
    assert(numelyx_core_numeric_softplus(DBL_MAX, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == DBL_MAX);
    assert(numelyx_core_numeric_softplus(-DBL_MAX, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.0);
    assert(numelyx_core_numeric_softplus(INFINITY, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 0.0);

    assert(numelyx_core_numeric_log_sigmoid(0.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert_close(value, -log(2.0), 1e-15);
    assert(numelyx_core_numeric_log_sigmoid(DBL_MAX, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.0);
    assert(numelyx_core_numeric_log_sigmoid(-DBL_MAX, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == -DBL_MAX);
    assert(numelyx_core_numeric_log_sigmoid(NAN, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == -DBL_MAX);
}

static void test_log_sum_exp_stability_and_validation(void)
{
    const double positive_values[] = {1000.0, 1000.0};
    const double negative_values[] = {-1000.0, -1001.0};
    const double reversed_values[] = {-1001.0, -1000.0};
    const double extreme_values[] = {DBL_MAX, -DBL_MAX};
    const double nonfinite_values[] = {1.0, INFINITY};
    const double singleton[] = {42.0};
    double value = 55.0;
    double reversed;

    assert(numelyx_core_numeric_log_sum_exp(
        positive_values, 2, &value
    ) == NUMERUS_NUMERIC_SUCCESS);
    assert_close(value, 1000.0 + log(2.0), 1e-12);

    assert(numelyx_core_numeric_log_sum_exp(
        negative_values, 2, &value
    ) == NUMERUS_NUMERIC_SUCCESS);
    assert_close(value, -1000.0 + log1p(exp(-1.0)), 1e-12);
    assert(numelyx_core_numeric_log_sum_exp(
        reversed_values, 2, &reversed
    ) == NUMERUS_NUMERIC_SUCCESS);
    assert_close(value, reversed, 1e-12);

    assert(numelyx_core_numeric_log_sum_exp(
        extreme_values, 2, &value
    ) == NUMERUS_NUMERIC_SUCCESS);
    assert(value == DBL_MAX);
    assert(numelyx_core_numeric_log_sum_exp(
        singleton, 1, &value
    ) == NUMERUS_NUMERIC_SUCCESS);
    assert(value == 42.0);

    value = 55.0;
    assert(numelyx_core_numeric_log_sum_exp(
        nonfinite_values, 2, &value
    ) == NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 55.0);
    assert(numelyx_core_numeric_log_sum_exp(
        NULL, 1, &value
    ) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
    assert(value == 55.0);
    assert(numelyx_core_numeric_log_sum_exp(
        singleton, 0, &value
    ) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
    assert(value == 55.0);
    assert(numelyx_core_numeric_log_sum_exp(
        singleton, 1, NULL
    ) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
}


static void test_log1mexp_stability_and_domain(void)
{
    const double near_zero = -1e-16;
    const double moderate = -1.0;
    double value = 91.0;

    assert(numelyx_core_numeric_log1mexp(0.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == -INFINITY);

    assert(numelyx_core_numeric_log1mexp(near_zero, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert_close(value, log(-expm1(near_zero)), 1e-15);

    assert(numelyx_core_numeric_log1mexp(moderate, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert_close(value, log1p(-exp(moderate)), 1e-15);

    value = 91.0;
    assert(numelyx_core_numeric_log1mexp(0.1, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 91.0);
    assert(numelyx_core_numeric_log1mexp(NAN, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 91.0);
    assert(numelyx_core_numeric_log1mexp(-1.0, NULL) ==
        NUMERUS_NUMERIC_INVALID_ARGUMENT);
}

static void test_log_diff_exp_stability_and_domain(void)
{
    const double a = 1000.0;
    const double b = 999.0;
    double value = 73.0;

    assert(numelyx_core_numeric_log_diff_exp(a, b, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert_close(value, a + log1p(-exp(b - a)), 1e-12);

    assert(numelyx_core_numeric_log_diff_exp(-2.0, -2.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == -INFINITY);

    assert(numelyx_core_numeric_log_diff_exp(DBL_MAX, -DBL_MAX, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == DBL_MAX);

    value = 73.0;
    assert(numelyx_core_numeric_log_diff_exp(1.0, 2.0, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 73.0);
    assert(numelyx_core_numeric_log_diff_exp(INFINITY, 0.0, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 73.0);
    assert(numelyx_core_numeric_log_diff_exp(1.0, 0.0, NULL) ==
        NUMERUS_NUMERIC_INVALID_ARGUMENT);
}

int main(void)
{
    test_sigmoid_stability();
    test_logit_domain_and_endpoints();
    test_softplus_and_log_sigmoid_extremes();
    test_log_sum_exp_stability_and_validation();
    test_log1mexp_stability_and_domain();
    test_log_diff_exp_stability_and_domain();
    puts("Stable scalar numerical tests passed.");
    return 0;
}
