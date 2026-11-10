/**
 * @file special_functions_test.c
 * @brief Reference and boundary tests for regularized incomplete functions.
 */
#include "../src/core/numelyx_core_special.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

static void close_to(double actual, double expected, double tolerance)
{
    if (!(fabs(actual - expected) <= tolerance)) {
        fprintf(stderr, "close_to failed: actual=%.17g expected=%.17g tolerance=%.17g\\n",
            actual, expected, tolerance);
        assert(fabs(actual - expected) <= tolerance);
    }
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


static void test_inverse_special_function_closed_forms(void)
{
    const double probabilities[] = {1e-8, 0.1, 0.5, 0.9, 1.0 - 1e-8};
    size_t i;
    for (i = 0; i < sizeof(probabilities) / sizeof(probabilities[0]); i++) {
        double value;
        double expected_p = -log1p(-probabilities[i]);
        double expected_q = -log(probabilities[i]);

        assert(numelyx_core_inverse_regularized_gamma_p(
            1.0, probabilities[i], &value) == NUMERUS_NUMERIC_SUCCESS);
        close_to(value, expected_p, 2e-13 * fmax(1.0, expected_p));

        assert(numelyx_core_inverse_regularized_gamma_q(
            1.0, probabilities[i], &value) == NUMERUS_NUMERIC_SUCCESS);
        close_to(value, expected_q, 2e-13 * fmax(1.0, expected_q));

        assert(numelyx_core_inverse_regularized_beta(
            1.0, 1.0, probabilities[i], &value) == NUMERUS_NUMERIC_SUCCESS);
        close_to(value, probabilities[i], 2e-14);
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


static void test_digamma_and_trigamma_reference_values(void)
{
    double value;
    /* psi(1) = -EulerGamma; psi(1/2) = -gamma - 2 log(2). */
    assert(numelyx_core_digamma(1.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, -0.5772156649015328606, 3e-14);
    assert(numelyx_core_digamma(0.5, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, -0.5772156649015328606 - 2.0 * log(2.0), 4e-14);
    assert(numelyx_core_digamma(8.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 2.0156414779556099965, 4e-14);

    /* psi_1(1) = pi^2/6; psi_1(1/2) = pi^2/2. */
    assert(numelyx_core_trigamma(1.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 1.6449340668482264365, 4e-14);
    assert(numelyx_core_trigamma(0.5, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 4.9348022005446793094, 8e-14);
    assert(numelyx_core_trigamma(8.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.1331370146940314251, 4e-14);
}

static void test_digamma_trigamma_validation_and_identities(void)
{
    double psi_x;
    double psi_next;
    double tri_x;
    double tri_next;
    double output = 23.0;

    assert(numelyx_core_digamma(2.25, &psi_x) == NUMERUS_NUMERIC_SUCCESS);
    assert(numelyx_core_digamma(3.25, &psi_next) == NUMERUS_NUMERIC_SUCCESS);
    close_to(psi_next - psi_x, 1.0 / 2.25, 2e-14);

    assert(numelyx_core_trigamma(2.25, &tri_x) == NUMERUS_NUMERIC_SUCCESS);
    assert(numelyx_core_trigamma(3.25, &tri_next) == NUMERUS_NUMERIC_SUCCESS);
    close_to(tri_x - tri_next, 1.0 / (2.25 * 2.25), 2e-14);

    assert(numelyx_core_digamma(0.0, &output) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(output == 23.0);
    assert(numelyx_core_trigamma(-1.0, &output) == NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(output == 23.0);
    assert(numelyx_core_digamma(INFINITY, &output) == NUMERUS_NUMERIC_NON_FINITE);
    assert(output == 23.0);
    assert(numelyx_core_trigamma(1.0, NULL) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
}


static void close_relative(double actual, double expected, double relative_tolerance)
{
    double scale = fabs(expected);
    double tolerance;
    if (scale < DBL_MIN) scale = DBL_MIN;
    tolerance = relative_tolerance * scale;
    assert(fabs(actual - expected) <= tolerance);
}

static void test_high_precision_gamma_reference_vectors(void)
{
    double value;
    assert(numelyx_core_regularized_gamma_p(0.1, 0.01, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.66262125995447979172, 3e-14);
    assert(numelyx_core_regularized_gamma_q(0.5, 0.2, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.52708925686553807367, 3e-14);
    assert(numelyx_core_regularized_gamma_p(100.0, 120.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.97213626010947933852, 4e-13);
    assert(numelyx_core_regularized_gamma_q(2.0, 50.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_relative(value, 9.8366242246159806934e-21, 2e-12);
}

static void test_high_precision_beta_reference_vectors(void)
{
    double value;
    assert(numelyx_core_regularized_beta(0.5, 0.5, 0.01, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.06376856085851984858, 3e-14);
    assert(numelyx_core_regularized_beta(2.5, 4.0, 0.3, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.35219758590676721388, 4e-14);
    assert(numelyx_core_regularized_beta(100.0, 120.0, 0.45, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 0.44780123014770762084, 4e-13);
}

static void test_special_function_monotonicity_and_digamma_derivative(void)
{
    double previous;
    double current;
    double psi_left;
    double psi_right;
    double trigamma;
    double x;
    size_t i;

    assert(numelyx_core_regularized_gamma_p(2.5, 0.1, &previous) ==
        NUMERUS_NUMERIC_SUCCESS);
    for (i = 2; i <= 40; i++) {
        x = (double) i / 10.0;
        assert(numelyx_core_regularized_gamma_p(2.5, x, &current) ==
            NUMERUS_NUMERIC_SUCCESS);
        assert(current >= previous);
        previous = current;
    }

    assert(numelyx_core_regularized_beta(2.5, 4.0, 0.01, &previous) ==
        NUMERUS_NUMERIC_SUCCESS);
    for (i = 2; i <= 99; i++) {
        x = (double) i / 100.0;
        assert(numelyx_core_regularized_beta(2.5, 4.0, x, &current) ==
            NUMERUS_NUMERIC_SUCCESS);
        assert(current >= previous);
        previous = current;
    }

    x = 2.375;
    assert(numelyx_core_digamma(x - 1e-5, &psi_left) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(numelyx_core_digamma(x + 1e-5, &psi_right) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(numelyx_core_trigamma(x, &trigamma) == NUMERUS_NUMERIC_SUCCESS);
    close_to((psi_right - psi_left) / 2e-5, trigamma, 2e-9);
}


static void test_erf_erfc_reference_values_and_symmetry(void)
{
    double erf_value;
    double erfc_value;
    double opposite;

    assert(numelyx_core_erf(0.0, &erf_value) == NUMERUS_NUMERIC_SUCCESS);
    assert(erf_value == 0.0);
    assert(numelyx_core_erf(1.0, &erf_value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(erf_value, 0.84270079294971486934, 3e-14);
    assert(numelyx_core_erf(-1.0, &erf_value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(erf_value, -0.84270079294971486934, 3e-14);

    assert(numelyx_core_erfc(1.0, &erfc_value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(erfc_value, 0.15729920705028513066, 3e-14);
    assert(numelyx_core_erfc(-1.0, &erfc_value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(erfc_value, 1.84270079294971486934, 3e-14);
    assert(numelyx_core_erfc(-1e-10, &erfc_value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(erfc_value, 1.0 + 1.1283791670955125739e-10, 3e-16);

    assert(numelyx_core_erfc(5.0, &erfc_value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(erfc_value, 1.5374597944280351e-12, 3e-25);
    assert(numelyx_core_erf(2.0, &erf_value) == NUMERUS_NUMERIC_SUCCESS);
    assert(numelyx_core_erf(-2.0, &opposite) == NUMERUS_NUMERIC_SUCCESS);
    close_to(erf_value, -opposite, 2e-15);
}

static void test_erf_erfc_validation_and_extremes(void)
{
    double value = 7.0;
    assert(numelyx_core_erf(NAN, &value) == NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 7.0);
    assert(numelyx_core_erfc(INFINITY, &value) == NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 7.0);
    assert(numelyx_core_erf(1.0, NULL) == NUMERUS_NUMERIC_INVALID_ARGUMENT);
    assert(numelyx_core_erf(DBL_MAX, &value) == NUMERUS_NUMERIC_SUCCESS);
    assert(value == 1.0);
    assert(numelyx_core_erfc(-DBL_MAX, &value) == NUMERUS_NUMERIC_SUCCESS);
    assert(value == 2.0);
    assert(numelyx_core_erfc(DBL_MAX, &value) == NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.0);
}


static void test_inverse_erf_erfc_reference_values(void)
{
    double value;
    assert(numelyx_core_inverse_erf(0.0, &value) == NUMERUS_NUMERIC_SUCCESS);
    assert(value == 0.0);
    assert(numelyx_core_inverse_erf(0.84270079294971486934, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 1.0, 3e-13);
    assert(numelyx_core_inverse_erf(-0.84270079294971486934, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, -1.0, 3e-13);
    assert(numelyx_core_inverse_erfc(0.15729920705028513066, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 1.0, 3e-13);
    assert(numelyx_core_inverse_erfc(1.84270079294971486934, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    close_to(value, -1.0, 3e-13);
    assert(numelyx_core_inverse_erfc(1.5, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, -0.47693627620446987338, 3e-13);
    assert(numelyx_core_inverse_erf(1e-10, &value) == NUMERUS_NUMERIC_SUCCESS);
    close_to(value, 8.8622692545275801365e-11, 2e-20);
}

static void test_normal_quantile_reference_values_and_round_trips(void)
{
    const double probabilities[] = {1e-10, 0.025, 0.5, 0.975, 1.0 - 1e-10};
    const double expected[] = {
        -6.361340902404056, -1.959963984540054,
        0.0, 1.959963984540054, 6.361340902404056
    };
    size_t i;
    for (i = 0; i < sizeof(probabilities) / sizeof(probabilities[0]); i++) {
        double quantile;
        double cdf;
        assert(numelyx_core_normal_quantile(probabilities[i], &quantile) ==
            NUMERUS_NUMERIC_SUCCESS);
        close_to(quantile, expected[i], 2e-10);
        assert(numelyx_core_erfc(fabs(quantile) / sqrt(2.0), &cdf) ==
            NUMERUS_NUMERIC_SUCCESS);
        cdf = quantile < 0.0 ? 0.5 * cdf : 1.0 - 0.5 * cdf;
        close_to(cdf, probabilities[i], 3e-13);
    }
}

static void test_inverse_error_function_validation(void)
{
    double value = 31.0;
    assert(numelyx_core_inverse_erf(1.1, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 31.0);
    assert(numelyx_core_inverse_erfc(-0.1, &value) ==
        NUMERUS_NUMERIC_DOMAIN_ERROR);
    assert(value == 31.0);
    assert(numelyx_core_normal_quantile(NAN, &value) ==
        NUMERUS_NUMERIC_NON_FINITE);
    assert(value == 31.0);
    assert(numelyx_core_normal_quantile(0.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == -INFINITY);
    assert(numelyx_core_normal_quantile(1.0, &value) ==
        NUMERUS_NUMERIC_SUCCESS);
    assert(value == INFINITY);
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
    test_inverse_special_function_closed_forms();
    test_digamma_and_trigamma_reference_values();
    test_digamma_trigamma_validation_and_identities();
    test_high_precision_gamma_reference_vectors();
    test_high_precision_beta_reference_vectors();
    test_special_function_monotonicity_and_digamma_derivative();
    test_erf_erfc_reference_values_and_symmetry();
    test_erf_erfc_validation_and_extremes();
    test_inverse_erf_erfc_reference_values();
    test_normal_quantile_reference_values_and_round_trips();
    test_inverse_error_function_validation();
    puts("Special-function tests passed.");
    return 0;
}
