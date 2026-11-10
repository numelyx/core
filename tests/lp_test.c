#include "../src/optimization/numelyx_core_lp.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static void close_to(double actual, double expected, double tolerance)
{
    assert(fabs(actual - expected) <= tolerance);
}

static void test_lp_reference_and_negative_rhs_phase_one(void)
{
    const double a[] = {
        1.0, 1.0,
        1.0, 0.0,
        0.0, 1.0
    };
    const double b[] = {4.0, 2.0, 3.0};
    const double c[] = {3.0, 2.0};
    const double phase_one_a[] = {-1.0, 1.0};
    const double phase_one_b[] = {-1.0, 3.0};
    const double phase_one_c[] = {1.0};
    numelyx_core_lp_options options;
    double solution[2] = {-9.0, -9.0};
    double objective = -9.0;
    size_t iterations = 0;

    assert(numelyx_core_lp_default_options(&options) == NUMELYX_LP_OPTIMAL);
    assert(numelyx_core_lp_solve_maximize(
        a, b, c, 3, 2, &options, solution, &objective, &iterations
    ) == NUMELYX_LP_OPTIMAL);
    close_to(solution[0], 2.0, 1e-8);
    close_to(solution[1], 2.0, 1e-8);
    close_to(objective, 10.0, 1e-8);
    assert(iterations > 0);

    solution[0] = -9.0;
    objective = -9.0;
    assert(numelyx_core_lp_solve_maximize(
        phase_one_a, phase_one_b, phase_one_c, 2, 1,
        &options, solution, &objective, NULL
    ) == NUMELYX_LP_OPTIMAL);
    close_to(solution[0], 3.0, 1e-8);
    close_to(objective, 3.0, 1e-8);
}

static void test_lp_infeasible_unbounded_and_zero_constraints(void)
{
    const double infeasible_a[] = {1.0, -1.0};
    const double infeasible_b[] = {0.0, -1.0};
    const double one_c[] = {1.0};
    const double unbounded_a[] = {1.0, -1.0};
    const double unbounded_b[] = {1.0};
    const double unbounded_c[] = {1.0, 1.0};
    const double negative_c[] = {-1.0, -2.0};
    numelyx_core_lp_options options;
    double solution[2] = {71.0, 72.0};
    double objective = 73.0;

    assert(numelyx_core_lp_default_options(&options) == NUMELYX_LP_OPTIMAL);
    assert(numelyx_core_lp_solve_maximize(
        infeasible_a, infeasible_b, one_c, 2, 1,
        &options, solution, &objective, NULL
    ) == NUMELYX_LP_INFEASIBLE);
    assert(solution[0] == 71.0 && objective == 73.0);

    assert(numelyx_core_lp_solve_maximize(
        unbounded_a, unbounded_b, unbounded_c, 1, 2,
        &options, solution, &objective, NULL
    ) == NUMELYX_LP_UNBOUNDED);
    assert(solution[0] == 71.0 && objective == 73.0);

    assert(numelyx_core_lp_solve_maximize(
        NULL, NULL, negative_c, 0, 2,
        &options, solution, &objective, NULL
    ) == NUMELYX_LP_OPTIMAL);
    assert(solution[0] == 0.0 && solution[1] == 0.0);
    assert(objective == 0.0);

    assert(numelyx_core_lp_solve_maximize(
        NULL, NULL, unbounded_c, 0, 2,
        &options, solution, &objective, NULL
    ) == NUMELYX_LP_UNBOUNDED);
}

static void test_lp_degenerate_and_iteration_limit(void)
{
    const double a[] = {
        1.0, 1.0,
        2.0, 2.0,
        1.0, 0.0,
        0.0, 1.0
    };
    const double b[] = {4.0, 8.0, 3.0, 3.0};
    const double c[] = {1.0, 1.0};
    numelyx_core_lp_options options;
    double solution[2] = {0.0, 0.0};
    double objective = 0.0;

    assert(numelyx_core_lp_default_options(&options) == NUMELYX_LP_OPTIMAL);
    options.max_iterations = 1;
    assert(numelyx_core_lp_solve_maximize(
        a, b, c, 4, 2, &options, solution, &objective, NULL
    ) == NUMELYX_LP_ITERATION_LIMIT);
}

static void test_lp_validation_and_output_preservation(void)
{
    const double a[] = {1.0};
    const double b[] = {1.0};
    const double c[] = {1.0};
    numelyx_core_lp_options options;
    double solution = 17.0;
    double objective = 18.0;

    assert(numelyx_core_lp_default_options(NULL) ==
        NUMELYX_LP_INVALID_ARGUMENT);
    assert(numelyx_core_lp_default_options(&options) == NUMELYX_LP_OPTIMAL);

    assert(numelyx_core_lp_solve_maximize(
        a, b, c, 1, 1, &options, NULL, &objective, NULL
    ) == NUMELYX_LP_INVALID_ARGUMENT);
    assert(solution == 17.0 && objective == 18.0);

    {
        const double nonfinite_c[] = {NAN};
        assert(numelyx_core_lp_solve_maximize(
            a, b, nonfinite_c, 1, 1, &options, &solution, &objective, NULL
        ) == NUMELYX_LP_NON_FINITE_INPUT);
        assert(solution == 17.0 && objective == 18.0);
    }

    {
        const double negative_a[] = {-INFINITY};
        assert(numelyx_core_lp_solve_maximize(
            negative_a, b, c, 1, 1, &options, &solution, &objective, NULL
        ) == NUMELYX_LP_NON_FINITE_INPUT);
        assert(solution == 17.0 && objective == 18.0);
    }

    options.feasibility_tolerance = 0.0;
    assert(numelyx_core_lp_solve_maximize(
        a, b, c, 1, 1, &options, &solution, &objective, NULL
    ) == NUMELYX_LP_INVALID_ARGUMENT);
    assert(solution == 17.0 && objective == 18.0);
}

int main(void)
{
    test_lp_reference_and_negative_rhs_phase_one();
    test_lp_infeasible_unbounded_and_zero_constraints();
    test_lp_degenerate_and_iteration_limit();
    test_lp_validation_and_output_preservation();
    puts("Linear programming tests passed.");
    return 0;
}
