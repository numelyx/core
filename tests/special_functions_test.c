/**
 * @file special_functions_test.c
 * @brief Reference and boundary tests for regularized incomplete functions.
 */
#include "../src/core/numelyx_core_special.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static void close_to(double actual, double expected, double tolerance)
{
    assert(fabs(actual - expected) <= tolerance);
}

static void test_regularized_gamma_reference_values(void)
{
    double value = -1.0;
    assert(numelyx_core_regularized_gamma_p(1.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 1.0 - exp(-1.0), 2e-14);

    assert(numelyx_core_regularized_gamma_q(1.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, exp(-1.0), 2e-14);

    assert(numelyx_core_regularized_gamma_p(0.5, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, erf(1.0), 2e-14);
}

static void test_regularized_gamma_tails_and_complementarity(void)
{
    double p;
    double q;

    assert(numelyx_core_regularized_gamma_p(1.0, 1e-8, &p) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(p, -expm1(-1e-8), 2e-15);

    assert(numelyx_core_regularized_gamma_q(1.0, 100.0, &q) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(q, exp(-100.0), 1e-57);

    assert(numelyx_core_regularized_gamma_p(2.5, 3.0, &p) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(numelyx_core_regularized_gamma_q(2.5, 3.0, &q) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(p + q, 1.0, 3e-15);
}

static void test_regularized_gamma_boundaries_and_validation(void)
{
    double value = 17.0;
    assert(numelyx_core_regularized_gamma_p(2.0, 0.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.0);
    assert(numelyx_core_regularized_gamma_q(2.0, 0.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == 1.0);

    value = 17.0;
    assert(numelyx_core_regularized_gamma_p(0.0, 1.0, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 17.0);
    assert(numelyx_core_regularized_gamma_q(1.0, -1.0, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 17.0);
    assert(numelyx_core_regularized_gamma_p(NAN, 1.0, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 17.0);
}

static void test_regularized_beta_reference_values(void)
{
    double value = -1.0;
    assert(numelyx_core_regularized_beta(2.0, 2.0, 0.5, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.5, 2e-14);

    assert(numelyx_core_regularized_beta(2.0, 3.0, 0.2, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.1808, 2e-14);
}

static void test_regularized_beta_symmetry_and_boundaries(void)
{
    double lower;
    double upper;
    double value = 11.0;

    assert(numelyx_core_regularized_beta(2.5, 4.0, 0.3, &lower) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(numelyx_core_regularized_beta(4.0, 2.5, 0.7, &upper) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(lower + upper, 1.0, 3e-14);

    assert(numelyx_core_regularized_beta(2.0, 3.0, 0.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.0);
    assert(numelyx_core_regularized_beta(2.0, 3.0, 1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == 1.0);

    value = 11.0;
    assert(numelyx_core_regularized_beta(2.0, 3.0, 1.1, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 11.0);
}

static void test_inverse_regularized_gamma_round_trips(void)
{
    const double probabilities[] = {1e-8, 0.1, 0.5, 0.9, 1.0 - 1e-8};
    size_t i;
    for (i = 0; i < sizeof(probabilities) / sizeof(probabilities[0]); i++) {
        double x;
        double p;
        double q;
        assert(numelyx_core_inverse_regularized_gamma_p(
            2.5, probabilities[i], &x) == NUMERUS_NUMERIC_SUCCESS);
        assert(numelyx_core_regularized_gamma_p(2.5, x, &p) ==
            NUMERUS_NUMERIC_SUCCESS);
        close_to(p, probabilities[i], 3e-13);

        assert(numelyx_core_inverse_regularized_gamma_q(
            2.5, probabilities[i], &x) == NUMERUS_NUMERIC_SUCCESS);
        assert(numelyx_core_regularized_gamma_q(2.5, x, &q) ==
            NUMERUS_NUMERIC_SUCCESS);
        close_to(q, probabilities[i], 3e-13);
    }
}

static void test_inverse_regularized_beta_round_trips(void)
{
    const double probabilities[] = {1e-8, 0.1, 0.5, 0.9, 1.0 - 1e-8};
    size_t i;
    for (i = 0; i < sizeof(probabilities) / sizeof(probabilities[0]); i++) {
        double x;
        double value;
        assert(numelyx_core_inverse_regularized_beta(
            2.5, 4.0, probabilities[i], &x) == NUMERUS_NUMERIC_SUCCESS);
        assert(numelyx_core_regularized_beta(2.5, 4.0, x, &value) ==
            NUMERUS_NUMERIC_SUCCESS);
        close_to(value, probabilities[i], 5e-13);
    }
}

static void test_inverse_special_function_validation(void)
{
    double value = 19.0;
    assert(numelyx_core_inverse_regularized_gamma_p(
        2.0, 0.0, &value) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 19.0);
    assert(numelyx_core_inverse_regularized_gamma_q(
        2.0, 1.0, &value) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 19.0);
    assert(numelyx_core_inverse_regularized_beta(
        2.0, 3.0, NAN, &value) == NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 19.0);
    assert(numelyx_core_inverse_regularized_gamma_p(
        2.0, 0.5, NULL) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
    assert(numelyx_core_inverse_regularized_gamma_q(
        2.0, 0.5, NULL) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
    assert(numelyx_core_inverse_regularized_beta(
        2.0, 3.0, 0.5, NULL) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
}

int main(void)
{
    test_regularized_gamma_reference_values();
    test_regularized_gamma_tails_and_complementarity();
    test_regularized_gamma_boundaries_and_validation();
    test_regularized_beta_reference_values();
    test_regularized_beta_symmetry_and_boundaries();
    test_inverse_regularized_gamma_round_trips();
    test_inverse_regularized_beta_round_trips();
    test_inverse_special_function_validation();
    puts("Special-function tests passed.");
    return 0;
}
