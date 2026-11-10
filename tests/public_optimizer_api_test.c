#include <math.h>

#include <numelyx/optimizer.h>

static bool objective(
    const double *parameters,
    size_t parameter_count,
    const void *context,
    double *value
)
{
    (void) context;
    if (parameters == NULL || parameter_count != 1 || value == NULL) {
        return false;
    }
    *value = (parameters[0] - 3.0) * (parameters[0] - 3.0);
    return true;
}

static bool gradient(
    const double *parameters,
    size_t parameter_count,
    const void *context,
    double *result
)
{
    (void) context;
    if (parameters == NULL || parameter_count != 1 || result == NULL) {
        return false;
    }
    result[0] = 2.0 * (parameters[0] - 3.0);
    return true;
}

int main(void)
{
    const double initial[] = {-1.0};
    double parameter = 0.0;
    double value = -1.0;
    numelyx_optimizer_options options;
    numelyx_optimizer_termination termination;
    numelyx_optimizer_result *result = NULL;

    if (numelyx_optimizer_default_options(&options) != NUMELYX_SUCCESS) {
        return 1;
    }
    if (numelyx_optimizer_minimize(
            objective,
            gradient,
            NULL,
            initial,
            1,
            &options,
            &result
        ) != NUMELYX_SUCCESS || result == NULL) {
        numelyx_optimizer_result_release(result);
        return 2;
    }

    if (numelyx_optimizer_result_parameter_count(result) != 1
        || numelyx_optimizer_result_get_parameter(result, 0, &parameter)
            != NUMELYX_SUCCESS
        || fabs(parameter - 3.0) > 1e-5
        || !numelyx_optimizer_result_objective_value_valid(result)
        || numelyx_optimizer_result_get_objective(result, &value)
            != NUMELYX_SUCCESS
        || value > 1e-8
        || numelyx_optimizer_result_get_termination(result, &termination)
            != NUMELYX_SUCCESS
        || termination != NUMELYX_OPTIMIZER_CONVERGED_GRADIENT) {
        numelyx_optimizer_result_release(result);
        return 3;
    }

    if (numelyx_optimizer_result_retain(result) != NUMELYX_SUCCESS) {
        numelyx_optimizer_result_release(result);
        return 4;
    }
    numelyx_optimizer_result_release(result);
    numelyx_optimizer_result_release(result);
    return 0;
}
