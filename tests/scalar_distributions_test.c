#include "../src/statistics/numelyx_core_distributions.h"
#include "../src/statistics/numelyx_core_rng.h"

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
    close_to(log_tail, -804.60844201375378817, 3e-12);
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

    value = 17.0;
    assert(numelyx_core_normal_pdf(0.0, 0.0, 0.0, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 17.0);
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


static void test_student_t_distribution_reference_values(void)
{
    double value;
    assert(numelyx_core_student_t_pdf(0.0, 1.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 1.0 / acos(-1.0), 3e-14);
    assert(numelyx_core_student_t_cdf(1.0, 1.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.75, 3e-14);
    assert(numelyx_core_student_t_survival(1.0, 1.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.25, 3e-14);
    assert(numelyx_core_student_t_quantile(0.75, 1.0, 0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 1.0, 3e-12);
    assert(numelyx_core_student_t_log_survival(
        10.0, 1.0, 0.0, 1.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    assert(isfinite(value) && value < -3.0);
}

static void test_chi_square_distribution_reference_values(void)
{
    double value;
    assert(numelyx_core_chi_square_pdf(2.0, 2.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.5 * exp(-1.0), 3e-14);
    assert(numelyx_core_chi_square_cdf(2.0, 2.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 1.0 - exp(-1.0), 3e-14);
    assert(numelyx_core_chi_square_survival(2.0, 2.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, exp(-1.0), 3e-14);
    assert(numelyx_core_chi_square_quantile(0.5, 2.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 2.0 * log(2.0), 3e-12);
    assert(numelyx_core_chi_square_log_survival(
        1000.0, 2.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, -500.0, 2e-12);
}

static void test_student_t_and_chi_square_validation(void)
{
    double value = 17.0;
    assert(numelyx_core_student_t_pdf(
        0.0, 0.0, 0.0, 1.0, &value) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 17.0);
    assert(numelyx_core_student_t_quantile(
        1.1, 2.0, 0.0, 1.0, &value) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 17.0);
    assert(numelyx_core_chi_square_pdf(
        1.0, 0.0, &value) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 17.0);
    assert(numelyx_core_chi_square_cdf(
        NAN, 2.0, &value) == NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 17.0);
}


static void test_normal_sampling_reproducibility_and_failure_atomicity(void)
{
    numelyx_core_rng *rng_a = NULL;
    numelyx_core_rng *rng_b = NULL;
    numelyx_core_rng *rng_c = NULL;
    double sample;
    double standard_normal;
    uint32_t next_a;
    uint32_t next_b;
    double untouched = 41.0;

    assert(numelyx_core_rng_create(UINT64_C(12345), UINT64_C(7), &rng_a) ==
        NUMERUS_RNG_SUCCESS);
    assert(numelyx_core_rng_clone(rng_a, &rng_b) == NUMERUS_RNG_SUCCESS);

    assert(numelyx_core_normal_sample(
        rng_a, 3.0, 2.0, &sample) == NUMERUS_NUMERIC_SUCCESS);
    assert(numelyx_core_rng_normal(rng_b, &standard_normal) ==
        NUMERUS_RNG_SUCCESS);
    close_to(sample, 3.0 + 2.0 * standard_normal, 2e-15);

    assert(numelyx_core_rng_next_u32(rng_a, &next_a) == NUMERUS_RNG_SUCCESS);
    assert(numelyx_core_rng_next_u32(rng_b, &next_b) == NUMERUS_RNG_SUCCESS);
    assert(next_a == next_b);

    assert(numelyx_core_normal_sample(
        rng_a, 0.0, 0.0, &untouched) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(untouched == 41.0);
    assert(numelyx_core_rng_next_u32(rng_a, &next_a) == NUMERUS_RNG_SUCCESS);
    assert(numelyx_core_rng_next_u32(rng_b, &next_b) == NUMERUS_RNG_SUCCESS);
    assert(next_a == next_b);

    assert(numelyx_core_rng_clone(rng_a, &rng_c) == NUMERUS_RNG_SUCCESS);
    assert(numelyx_core_rng_normal(rng_c, &standard_normal) ==
        NUMERUS_RNG_SUCCESS);
    untouched = 41.0;
    assert(numelyx_core_normal_sample(
        rng_a,
        copysign(DBL_MAX, standard_normal),
        DBL_MAX,
        &untouched
    ) == NUMERUS_NUMERIC_NUMERICAL_FAILURE);
    assert(untouched == 41.0);
    assert(numelyx_core_rng_next_u32(rng_a, &next_a) == NUMERUS_RNG_SUCCESS);
    assert(numelyx_core_rng_next_u32(rng_c, &next_b) == NUMERUS_RNG_SUCCESS);
    assert(next_a == next_b);

    numelyx_core_rng_destroy(rng_c);
    numelyx_core_rng_destroy(rng_b);
    numelyx_core_rng_destroy(rng_a);
}

int main(void)
{
    test_normal_distribution_reference_values();
    test_normal_distribution_log_tails();
    test_normal_distribution_quantile_and_validation();
    test_student_t_distribution_reference_values();
    test_chi_square_distribution_reference_values();
    test_student_t_and_chi_square_validation();
    test_normal_sampling_reproducibility_and_failure_atomicity();
    puts("Scalar distribution tests passed.");
    return 0;
}
