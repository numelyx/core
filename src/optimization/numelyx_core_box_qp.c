#include "numelyx_core_box_qp.h"
#include "numelyx_core_size.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

static double clamp_box(double x, double lower, double upper)
{
    return x < lower ? lower : (x > upper ? upper : x);
}

static bool qp_objective(const double *h, const double *g, const double *x,
    size_t n, double *value)
{
    size_t i, j;
    double quadratic = 0.0, linear = 0.0;
    for (i = 0; i < n; i++) {
        double row = 0.0;
        for (j = 0; j < n; j++) {
            row += h[i*n+j] * x[j];
            if (!isfinite(row)) return false;
        }
        quadratic += x[i] * row;
        linear += g[i] * x[i];
        if (!isfinite(quadratic) || !isfinite(linear)) return false;
    }
    *value = 0.5 * quadratic + linear;
    return isfinite(*value);
}

static bool qp_gradient(const double *h, const double *g, const double *x,
    size_t n, double *gradient)
{
    size_t i, j;
    for (i = 0; i < n; i++) {
        double value = g[i];
        for (j = 0; j < n; j++) {
            value += h[i*n+j] * x[j];
            if (!isfinite(value)) return false;
        }
        gradient[i] = value;
    }
    return true;
}

static bool qp_validate_psd(const double *h, size_t n, double tolerance,
    double *factor, double *lipschitz)
{
    size_t i, j, k;
    double scale = 0.0, max_row_sum = 0.0;
    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) {
            double value = h[i*n+j];
            if (!isfinite(value)) return false;
            if (fabs(value) > scale) scale = fabs(value);
            row_sum += fabs(value);
            if (!isfinite(row_sum)) return false;
        }
        if (row_sum > max_row_sum) max_row_sum = row_sum;
    }
    for (i = 0; i < n; i++) {
        for (j = i+1; j < n; j++) {
            if (fabs(h[i*n+j]-h[j*n+i]) > tolerance*scale) return false;
        }
    }
    memset(factor, 0, n*n*sizeof(*factor));
    for (i = 0; i < n; i++) {
        double diagonal = h[i*n+i];
        for (k = 0; k < i; k++) {
            double v = factor[i*n+k];
            diagonal -= v*v;
        }
        if (!isfinite(diagonal) || diagonal < -tolerance*scale) return false;
        if (diagonal <= tolerance*scale) {
            factor[i*n+i] = 0.0;
            for (j = i+1; j < n; j++) {
                double residual = h[j*n+i];
                for (k = 0; k < i; k++)
                    residual -= factor[j*n+k]*factor[i*n+k];
                if (!isfinite(residual) || fabs(residual) > tolerance*scale)
                    return false;
                factor[j*n+i] = 0.0;
            }
        } else {
            factor[i*n+i] = sqrt(diagonal);
            for (j = i+1; j < n; j++) {
                double residual = h[j*n+i];
                for (k = 0; k < i; k++)
                    residual -= factor[j*n+k]*factor[i*n+k];
                factor[j*n+i] = residual/factor[i*n+i];
                if (!isfinite(factor[j*n+i])) return false;
            }
        }
    }
    *lipschitz = max_row_sum;
    return true;
}

numelyx_core_qp_status numelyx_core_box_qp_default_options(
    numelyx_core_box_qp_options *options)
{
    if (options == NULL) return NUMELYX_QP_INVALID_ARGUMENT;
    options->max_iterations = 20000;
    options->projected_gradient_tolerance = 1e-8;
    options->positive_semidefinite_tolerance = 1e-10;
    return NUMELYX_QP_SUCCESS;
}

numelyx_core_qp_status numelyx_core_box_qp_solve(
    const double *h, const double *g, const double *lower,
    const double *upper, const double *initial, size_t n,
    const numelyx_core_box_qp_options *options, double *solution,
    double *objective_value, size_t *iterations)
{
    size_t matrix_elements, matrix_bytes, vector_bytes, i, iteration;
    double *factor=NULL, *x=NULL, *y=NULL, *candidate=NULL;
    double *gradient=NULL, *candidate_gradient=NULL;
    double lipschitz=0.0, current_objective, momentum=1.0;
    numelyx_core_qp_status status=NUMELYX_QP_NUMERICAL_FAILURE;

    if (h==NULL || g==NULL || lower==NULL || upper==NULL || initial==NULL ||
        solution==NULL || objective_value==NULL || options==NULL || n==0 ||
        options->max_iterations==0) return NUMELYX_QP_INVALID_ARGUMENT;
    if (!isfinite(options->projected_gradient_tolerance) ||
        !isfinite(options->positive_semidefinite_tolerance))
        return NUMELYX_QP_NON_FINITE_INPUT;
    if (options->projected_gradient_tolerance<=0.0 ||
        options->positive_semidefinite_tolerance<=0.0)
        return NUMELYX_QP_INVALID_ARGUMENT;
    if (!numelyx_core_size_multiply(n,n,&matrix_elements) ||
        !numelyx_core_size_multiply(matrix_elements,sizeof(*factor),&matrix_bytes) ||
        !numelyx_core_size_multiply(n,sizeof(*x),&vector_bytes))
        return NUMELYX_QP_SIZE_OVERFLOW;

    for (i=0;i<matrix_elements;i++) {
        if (!isfinite(h[i])) return NUMELYX_QP_NON_FINITE_INPUT;
    }
    for (i=0;i<n;i++) {
        if (!isfinite(g[i]) || !isfinite(initial[i]) ||
            isnan(lower[i]) || isnan(upper[i]))
            return NUMELYX_QP_NON_FINITE_INPUT;
        if (lower[i]==INFINITY || upper[i]==-INFINITY ||
            lower[i]>upper[i])
            return NUMELYX_QP_INVALID_ARGUMENT;
    }
    factor=(double*)calloc(1,matrix_bytes);
    x=(double*)malloc(vector_bytes);
    y=(double*)malloc(vector_bytes);
    candidate=(double*)malloc(vector_bytes);
    gradient=(double*)malloc(vector_bytes);
    candidate_gradient=(double*)malloc(vector_bytes);
    if (!factor || !x || !y || !candidate || !gradient || !candidate_gradient) {
        status=NUMELYX_QP_OUT_OF_MEMORY;
        goto cleanup;
    }
    if (!qp_validate_psd(h,n,options->positive_semidefinite_tolerance,
        factor,&lipschitz)) {
        status=NUMELYX_QP_NOT_POSITIVE_SEMIDEFINITE;
        goto cleanup;
    }
    for (i=0;i<n;i++) x[i]=clamp_box(initial[i],lower[i],upper[i]);

    if (lipschitz==0.0) {
        for (i=0;i<n;i++) {
            if (g[i]>0.0) {
                if (!isfinite(lower[i])) { status=NUMELYX_QP_UNBOUNDED; goto cleanup; }
                x[i]=lower[i];
            } else if (g[i]<0.0) {
                if (!isfinite(upper[i])) { status=NUMELYX_QP_UNBOUNDED; goto cleanup; }
                x[i]=upper[i];
            }
        }
        if (!qp_objective(h,g,x,n,&current_objective)) goto cleanup;
        memcpy(solution,x,vector_bytes);
        *objective_value=current_objective;
        if (iterations) *iterations=0;
        status=NUMELYX_QP_CONVERGED;
        goto cleanup;
    }

    memcpy(y,x,vector_bytes);
    if (!qp_objective(h,g,x,n,&current_objective)) goto cleanup;

    for (iteration=1;iteration<=options->max_iterations;iteration++) {
        double trial_objective, maximum_mapping=0.0, scale=1.0;
        double next_momentum, extrapolation;
        bool restarted=false;
        if (!qp_gradient(h,g,y,n,gradient)) goto cleanup;
        for (i=0;i<n;i++) {
            double point=y[i]-gradient[i]/lipschitz;
            if (!isfinite(point)) goto cleanup;
            candidate[i]=clamp_box(point,lower[i],upper[i]);
        }
        if (!qp_objective(h,g,candidate,n,&trial_objective)) goto cleanup;
        if (trial_objective > current_objective +
            32.0*DBL_EPSILON*fmax(1.0,fabs(current_objective))) {
            restarted=true;
            momentum=1.0;
            memcpy(y,x,vector_bytes);
            if (!qp_gradient(h,g,y,n,gradient)) goto cleanup;
            for (i=0;i<n;i++) {
                double point=x[i]-gradient[i]/lipschitz;
                if (!isfinite(point)) goto cleanup;
                candidate[i]=clamp_box(point,lower[i],upper[i]);
            }
            if (!qp_objective(h,g,candidate,n,&trial_objective) ||
                trial_objective > current_objective +
                64.0*DBL_EPSILON*fmax(1.0,fabs(current_objective))) goto cleanup;
        }
        if (!qp_gradient(h,g,candidate,n,candidate_gradient)) goto cleanup;
        for (i=0;i<n;i++) {
            double projected=clamp_box(
                candidate[i]-candidate_gradient[i]/lipschitz,lower[i],upper[i]);
            double mapping=lipschitz*(candidate[i]-projected);
            if (!isfinite(mapping)) goto cleanup;
            if (fabs(mapping)>maximum_mapping) maximum_mapping=fabs(mapping);
            if (fabs(candidate[i])>scale) scale=fabs(candidate[i]);
        }
        if (maximum_mapping <= options->projected_gradient_tolerance*scale) {
            memcpy(solution,candidate,vector_bytes);
            *objective_value=trial_objective;
            if (iterations) *iterations=iteration;
            status=NUMELYX_QP_CONVERGED;
            goto cleanup;
        }
        next_momentum=0.5*(1.0+sqrt(1.0+4.0*momentum*momentum));
        extrapolation=(momentum-1.0)/next_momentum;
        for (i=0;i<n;i++) {
            y[i]=candidate[i]+extrapolation*(candidate[i]-x[i]);
            if (!isfinite(y[i])) goto cleanup;
        }
        memcpy(x,candidate,vector_bytes);
        current_objective=trial_objective;
        momentum=next_momentum;
        if (restarted) {
            momentum=1.0;
            memcpy(y,x,vector_bytes);
        }
    }
    memcpy(solution,x,vector_bytes);
    *objective_value=current_objective;
    if (iterations) *iterations=options->max_iterations;
    status=NUMELYX_QP_ITERATION_LIMIT;

cleanup:
    free(factor); free(x); free(y); free(candidate);
    free(gradient); free(candidate_gradient);
    return status;
}
