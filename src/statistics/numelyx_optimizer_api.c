#include <numelyx/optimizer.h>

#include "numelyx_core_optimizer.h"

static numelyx_status map_optimizer_status(numelyx_core_optimizer_status status)
{
    switch (status) {
        case NUMERUS_OPTIMIZER_SUCCESS:
            return NUMELYX_SUCCESS;
        case NUMERUS_OPTIMIZER_INVALID_ARGUMENT:
            return NUMELYX_INVALID_ARGUMENT;
        case NUMERUS_OPTIMIZER_NON_FINITE_INPUT:
            return NUMELYX_NON_FINITE;
        case NUMERUS_OPTIMIZER_SIZE_OVERFLOW:
            return NUMELYX_SIZE_OVERFLOW;
        case NUMERUS_OPTIMIZER_OUT_OF_MEMORY:
            return NUMELYX_OUT_OF_MEMORY;
        case NUMERUS_OPTIMIZER_OUT_OF_BOUNDS:
            return NUMELYX_OUT_OF_BOUNDS;
        case NUMERUS_OPTIMIZER_NO_OBJECTIVE_VALUE:
        default:
            return NUMELYX_NUMERICAL_ERROR;
    }
}

static numelyx_optimizer_termination map_termination(
    numelyx_core_optimizer_termination termination
)
{
    switch (termination) {
        case NUMERUS_OPTIMIZER_CONVERGED_GRADIENT:
            return NUMELYX_OPTIMIZER_CONVERGED_GRADIENT;
        case NUMERUS_OPTIMIZER_ITERATION_LIMIT:
            return NUMELYX_OPTIMIZER_ITERATION_LIMIT;
        case NUMERUS_OPTIMIZER_STEP_TOO_SMALL:
            return NUMELYX_OPTIMIZER_STEP_TOO_SMALL;
        case NUMERUS_OPTIMIZER_LINE_SEARCH_FAILED:
            return NUMELYX_OPTIMIZER_LINE_SEARCH_FAILED;
        case NUMERUS_OPTIMIZER_OBJECTIVE_CALLBACK_FAILED:
            return NUMELYX_OPTIMIZER_OBJECTIVE_CALLBACK_FAILED;
        case NUMERUS_OPTIMIZER_GRADIENT_CALLBACK_FAILED:
            return NUMELYX_OPTIMIZER_GRADIENT_CALLBACK_FAILED;
        case NUMERUS_OPTIMIZER_GRADIENT_EVALUATION_FAILED:
            return NUMELYX_OPTIMIZER_GRADIENT_EVALUATION_FAILED;
        case NUMERUS_OPTIMIZER_NON_FINITE_OBJECTIVE:
            return NUMELYX_OPTIMIZER_NON_FINITE_OBJECTIVE;
        case NUMERUS_OPTIMIZER_NON_FINITE_GRADIENT:
            return NUMELYX_OPTIMIZER_NON_FINITE_GRADIENT;
        case NUMERUS_OPTIMIZER_NUMERICAL_BREAKDOWN:
        default:
            return NUMELYX_OPTIMIZER_NUMERICAL_BREAKDOWN;
    }
}

static void copy_options_to_core(
    const numelyx_optimizer_options *source,
    numelyx_core_optimizer_options *destination
)
{
    destination->max_iterations = source->max_iterations;
    destination->max_line_search_iterations = source->max_line_search_iterations;
    destination->gradient_tolerance = source->gradient_tolerance;
    destination->minimum_step = source->minimum_step;
    destination->step_tolerance = source->step_tolerance;
    destination->initial_step = source->initial_step;
    destination->finite_difference_relative_step =
        source->finite_difference_relative_step;
}

numelyx_status numelyx_optimizer_default_options(
    numelyx_optimizer_options *options
)
{
    numelyx_core_optimizer_options core_options;
    numelyx_core_optimizer_status status;

    if (options == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }

    status = numelyx_core_optimizer_default_options(&core_options);
    if (status != NUMERUS_OPTIMIZER_SUCCESS) {
        return map_optimizer_status(status);
    }

    options->struct_size = (uint32_t) sizeof(*options);
    options->version = NUMELYX_OPTIMIZER_OPTIONS_VERSION_1;
    options->max_iterations = core_options.max_iterations;
    options->max_line_search_iterations = core_options.max_line_search_iterations;
    options->gradient_tolerance = core_options.gradient_tolerance;
    options->minimum_step = core_options.minimum_step;
    options->step_tolerance = core_options.step_tolerance;
    options->initial_step = core_options.initial_step;
    options->finite_difference_relative_step =
        core_options.finite_difference_relative_step;

    return NUMELYX_SUCCESS;
}

numelyx_status numelyx_optimizer_minimize(
    numelyx_objective_fn objective,
    numelyx_gradient_fn gradient,
    const void *context,
    const double *initial_parameters,
    size_t parameter_count,
    const numelyx_optimizer_options *options,
    numelyx_optimizer_result **result
)
{
    numelyx_core_optimizer_options core_options;
    const numelyx_core_optimizer_options *core_options_pointer = NULL;

    if (result == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    *result = NULL;

    if (options == NULL
        || options->struct_size < (uint32_t) sizeof(*options)
        || options->version != NUMELYX_OPTIMIZER_OPTIONS_VERSION_1) {
        return NUMELYX_INVALID_ARGUMENT;
    }
    copy_options_to_core(options, &core_options);
    core_options_pointer = &core_options;

    return map_optimizer_status(numelyx_core_optimizer_minimize(
        objective,
        gradient,
        context,
        initial_parameters,
        parameter_count,
        core_options_pointer,
        result
    ));
}

numelyx_status numelyx_optimizer_result_retain(
    numelyx_optimizer_result *result
)
{
    return map_optimizer_status(numelyx_core_optimizer_result_retain(result));
}

void numelyx_optimizer_result_release(numelyx_optimizer_result *result)
{
    numelyx_core_optimizer_result_destroy(result);
}

size_t numelyx_optimizer_result_parameter_count(
    const numelyx_optimizer_result *result
)
{
    return numelyx_core_optimizer_result_parameter_count(result);
}

numelyx_status numelyx_optimizer_result_get_parameter(
    const numelyx_optimizer_result *result,
    size_t index,
    double *value
)
{
    return map_optimizer_status(
        numelyx_core_optimizer_result_get_parameter(result, index, value)
    );
}

bool numelyx_optimizer_result_objective_value_valid(
    const numelyx_optimizer_result *result
)
{
    return numelyx_core_optimizer_result_objective_value_valid(result);
}

numelyx_status numelyx_optimizer_result_get_objective(
    const numelyx_optimizer_result *result,
    double *value
)
{
    return map_optimizer_status(
        numelyx_core_optimizer_result_get_objective(result, value)
    );
}

numelyx_status numelyx_optimizer_result_get_termination(
    const numelyx_optimizer_result *result,
    numelyx_optimizer_termination *termination
)
{
    numelyx_core_optimizer_termination core_termination;
    numelyx_core_optimizer_status status;

    if (termination == NULL) {
        return NUMELYX_INVALID_ARGUMENT;
    }

    status = numelyx_core_optimizer_result_get_termination(
        result,
        &core_termination
    );
    if (status != NUMERUS_OPTIMIZER_SUCCESS) {
        return map_optimizer_status(status);
    }

    *termination = map_termination(core_termination);
    return NUMELYX_SUCCESS;
}

size_t numelyx_optimizer_result_iterations(
    const numelyx_optimizer_result *result
)
{
    return numelyx_core_optimizer_result_iterations(result);
}

size_t numelyx_optimizer_result_objective_evaluations(
    const numelyx_optimizer_result *result
)
{
    return numelyx_core_optimizer_result_objective_evaluations(result);
}

size_t numelyx_optimizer_result_gradient_evaluations(
    const numelyx_optimizer_result *result
)
{
    return numelyx_core_optimizer_result_gradient_evaluations(result);
}
