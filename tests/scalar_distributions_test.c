#include "../src/core/numelyx_core_distributions.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

static void close_to(double actual, double expected, double tolerance)
{
    assert(fabs(actual - expected) <= tolerance);
}

static void test_normal_distribution_reference_values(void)
{
    double value;
    assert(numelyx_core_normal_log_pdf(0.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, -0.91893853320467274178, 2e-15);
    assert(numelyx_core_normal_pdf(0.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.39894228040143267794, 2e-15);

    assert(numelyx_core_normal_cdf(0.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.5, 2e-15);
    assert(numelyx_core_normal_survival(0.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.5, 2e-15);

    assert(numelyx_core_normal_cdf(1.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.84134474606854294859, 3e-14);
    assert(numelyx_core_normal_survival(1.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.15865525393145705141, 3e-14);
}

static void test_normal_distribution_log_tails(void)
{
    double log_tail;
    double tail;
    assert(numelyx_core_normal_log_cdf(-10.0, 0.0, 1.0, &log_tail) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(log_tail, -53.23128515051247058, 2e-12);
    assert(numelyx_core_normal_log_survival(10.0, 0.0, 1.0, &log_tail) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(log_tail, -53.23128515051247058, 2e-12);
    assert(numelyx_core_normal_survival(10.0, 0.0, 1.0, &tail) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(tail, 7.619853024160526e-24, 2e-36);

    assert(numelyx_core_normal_log_cdf(-40.0, 0.0, 1.0, &log_tail) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(isfinite(log_tail) && log_tail < -800.0);
}

static void test_normal_distribution_quantile_and_validation(void)
{
    double value = 19.0;
    assert(numelyx_core_normal_distribution_quantile(
        0.975, 0.0, 1.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 1.959963984540054, 2e-12);
    assert(numelyx_core_normal_distribution_quantile(
        0.025, 3.0, 2.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 3.0 - 2.0 * 1.959963984540054, 4e-12);

    assert(numelyx_core_normal_pdf(0.0, 0.0, 0.0, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 3.0 - 2.0 * 1.959963984540054);
    assert(numelyx_core_normal_cdf(NAN, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(numelyx_core_normal_distribution_quantile(
        1.1, 0.0, 1.0, &value) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(numelyx_core_normal_distribution_quantile(
        0.5, 0.0, 1.0, NULL) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
    assert(numelyx_core_normal_log_pdf(DBL_MAX, -DBL_MAX, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == -INFINITY);
}

int main(void)
{
    test_normal_distribution_reference_values();
    test_normal_distribution_log_tails();
    test_normal_distribution_quantile_and_validation();
    puts("Scalar distribution tests passed.");
    return 0;
}
