#include <float.h>
#include <math.h>
#include <stddef.h>

#include <numelyx/matrix.h>

static int close_to(double actual, double expected)
{
    return fabs(actual - expected) < 1e-8;
}

int main(void)
{
    numelyx_matrix *a = NULL, *b = NULL, *q = NULL, *r = NULL;
    numelyx_matrix *u = NULL, *s = NULL, *vt = NULL, *pinv = NULL;
    numelyx_matrix *lower = NULL, *solution = NULL, *triangular = NULL;
    numelyx_matrix *rhs2 = NULL, *tri_solution = NULL;
    numelyx_matrix *basis = NULL, *extra = NULL;
    size_t rank = 0, subspace_dimension = 99;
    numelyx_matrix_solution_kind classification;
    int sign = 0;
    double scalar = 0.0, log_abs = 0.0;
    const double diagonal[] = {2.0, 0.0, 0.0, 3.0};
    const double design[] = {1.0, 0.0, 1.0, 1.0, 1.0, 2.0};
    const double response[] = {1.0, 3.0, 5.0};
    const double weights[] = {1.0, 2.0, 1.0};
    const double spd[] = {4.0, 2.0, 2.0, 3.0};
    const double tri[] = {2.0, 0.0, 1.0, 3.0};
    const double rhs_values[] = {2.0, 4.0};

    if (numelyx_matrix_create_dense(2, 2, diagonal, &a) != NUMELYX_SUCCESS)
        goto fail;
    if (numelyx_matrix_log_determinant(a, &sign, &log_abs) != NUMELYX_SUCCESS
        || sign != 1 || !close_to(log_abs, log(6.0)))
        goto fail;
    if (numelyx_matrix_rank(a, &rank) != NUMELYX_SUCCESS || rank != 2)
        goto fail;
    if (numelyx_matrix_column_space_basis(
            a, &basis, &subspace_dimension) != NUMELYX_SUCCESS
        || basis == NULL || subspace_dimension != 2)
        goto fail;
    numelyx_matrix_release(basis); basis = NULL;
    if (numelyx_matrix_row_space_basis(
            a, &basis, &subspace_dimension) != NUMELYX_SUCCESS
        || basis == NULL || subspace_dimension != 2)
        goto fail;
    numelyx_matrix_release(basis); basis = NULL;
    subspace_dimension = 99;
    if (numelyx_matrix_null_space(
            a, &basis, &subspace_dimension) != NUMELYX_SUCCESS
        || basis != NULL || subspace_dimension != 0)
        goto fail;
    {
        int has_nan = 1, has_infinity = 1, is_finite = 0;
        if (numelyx_matrix_has_nan(a, &has_nan) != NUMELYX_SUCCESS
            || has_nan != 0
            || numelyx_matrix_has_infinity(a, &has_infinity)
                != NUMELYX_SUCCESS || has_infinity != 0
            || numelyx_matrix_is_finite(a, &is_finite)
                != NUMELYX_SUCCESS || is_finite != 1)
            goto fail;
    }
    if (numelyx_matrix_trace(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 5.0))
        goto fail;
    if (numelyx_matrix_norm_frobenius(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, sqrt(13.0)))
        goto fail;
    if (numelyx_matrix_norm_one(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 3.0))
        goto fail;
    if (numelyx_matrix_norm_infinity(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 3.0))
        goto fail;

    if (numelyx_matrix_qr_decompose(a, &q, &q) != NUMELYX_INVALID_ARGUMENT
        || q != NULL)
        goto fail;
    if (numelyx_matrix_qr_decompose(a, &q, &r) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(q) != 2 || numelyx_matrix_columns(q) != 2
        || numelyx_matrix_rows(r) != 2 || numelyx_matrix_columns(r) != 2)
        goto fail;
    if (numelyx_matrix_svd(a, &u, &u, &vt) != NUMELYX_INVALID_ARGUMENT
        || u != NULL || vt != NULL)
        goto fail;
    if (numelyx_matrix_svd(a, &u, &s, &vt) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(u) != 2
        || numelyx_matrix_rows(s) != 2
        || numelyx_matrix_columns(s) != 1
        || numelyx_matrix_columns(vt) != 2)
        goto fail;
    if (numelyx_matrix_pseudoinverse(a, &pinv) != NUMELYX_SUCCESS
        || numelyx_matrix_get(pinv, 0, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 0.5)
        || numelyx_matrix_get(pinv, 1, 1, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 1.0 / 3.0))
        goto fail;

    if (numelyx_matrix_create_dense(2, 2, spd, &b) != NUMELYX_SUCCESS
        || numelyx_matrix_cholesky(b, &lower) != NUMELYX_SUCCESS
        || numelyx_matrix_get(lower, 0, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 2.0))
        goto fail;

    if (numelyx_matrix_symmetric_square_root(
            b, &extra) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(extra) != 2)
        goto fail;
    numelyx_matrix_release(extra); extra = NULL;
    if (numelyx_matrix_symmetric_logarithm(
            b, &extra) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(extra) != 2)
        goto fail;
    numelyx_matrix_release(extra); extra = NULL;
    if (numelyx_matrix_symmetric_sine(
            b, &extra) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(extra) != 2)
        goto fail;
    numelyx_matrix_release(extra); extra = NULL;
    if (numelyx_matrix_symmetric_cosine(
            b, &extra) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(extra) != 2)
        goto fail;
    numelyx_matrix_release(extra); extra = NULL;
    if (numelyx_matrix_row_echelon_form(
            a, &extra) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(extra) != 2)
        goto fail;
    numelyx_matrix_release(extra); extra = NULL;
    if (numelyx_matrix_reduced_row_echelon_form(
            a, &extra) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(extra) != 2)
        goto fail;
    numelyx_matrix_release(extra); extra = NULL;
    {
        int skew = 1, orthogonal = 1;
        if (numelyx_matrix_symmetric_spectral_radius(a, &scalar)
                != NUMELYX_SUCCESS || !close_to(scalar, 3.0)
            || numelyx_matrix_is_skew_symmetric(a, &skew)
                != NUMELYX_SUCCESS || skew != 0
            || numelyx_matrix_is_orthogonal(a, &orthogonal)
                != NUMELYX_SUCCESS || orthogonal != 0)
            goto fail;
    }

    numelyx_matrix_release(u); u = NULL;
    numelyx_matrix_release(vt); vt = NULL;
    if (numelyx_matrix_symmetric_eigen(a, &u, &u)
            != NUMELYX_INVALID_ARGUMENT || u != NULL)
        goto fail;
    if (numelyx_matrix_symmetric_eigen(
            a, &u, &vt) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(u) != 2
        || numelyx_matrix_columns(u) != 1
        || numelyx_matrix_rows(vt) != 2
        || numelyx_matrix_columns(vt) != 2
        || numelyx_matrix_get(u, 0, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 3.0)
        || numelyx_matrix_get(u, 1, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 2.0))
        goto fail;
    numelyx_matrix_release(u); u = NULL;
    numelyx_matrix_release(vt); vt = NULL;
    if (numelyx_matrix_ldlt_decompose(a, &u, &u)
            != NUMELYX_INVALID_ARGUMENT || u != NULL)
        goto fail;
    if (numelyx_matrix_ldlt_decompose(
            a, &u, &vt) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(u) != 2
        || numelyx_matrix_rows(vt) != 2
        || numelyx_matrix_get(u, 0, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 1.0)
        || numelyx_matrix_get(u, 0, 1, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 0.0)
        || numelyx_matrix_get(u, 1, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 0.0)
        || numelyx_matrix_get(u, 1, 1, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 1.0)
        || numelyx_matrix_get(vt, 0, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 2.0)
        || numelyx_matrix_get(vt, 0, 1, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 0.0)
        || numelyx_matrix_get(vt, 1, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 0.0)
        || numelyx_matrix_get(vt, 1, 1, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 3.0))
        goto fail;
    numelyx_matrix_release(u); u = NULL;
    numelyx_matrix_release(vt); vt = NULL;
    if (numelyx_matrix_spectral_norm(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 3.0)
        || numelyx_matrix_condition_estimate_one(a, &scalar)
            != NUMELYX_SUCCESS || !close_to(scalar, 1.5))
        goto fail;
    {
        const double pivoted_values[] = {0.0, 1.0, 2.0, 3.0};
        const double singular_values[] = {2.0, 0.0, 0.0, 0.0};
        numelyx_matrix *pivoted = NULL;
        numelyx_matrix *singular = NULL;

        if (numelyx_matrix_create_dense(
                2, 2, pivoted_values, &pivoted) != NUMELYX_SUCCESS
            || numelyx_matrix_condition_estimate_one(pivoted, &scalar)
                != NUMELYX_SUCCESS || !close_to(scalar, 10.0)
            || numelyx_matrix_create_dense(
                2, 2, singular_values, &singular) != NUMELYX_SUCCESS
            || numelyx_matrix_condition_estimate_one(singular, &scalar)
                != NUMELYX_SUCCESS || !isinf(scalar)) {
            numelyx_matrix_release(singular);
            numelyx_matrix_release(pivoted);
            goto fail;
        }

        numelyx_matrix_release(singular);
        numelyx_matrix_release(pivoted);
    }
    {
        const double tiny_value[] = {DBL_MIN};
        const double huge_value[] = {DBL_MAX};
        numelyx_matrix *tiny = NULL;
        numelyx_matrix *huge = NULL;

        if (numelyx_matrix_create_dense(1, 1, tiny_value, &tiny)
                != NUMELYX_SUCCESS
            || numelyx_matrix_condition_estimate_one(tiny, &scalar)
                != NUMELYX_SUCCESS || !close_to(scalar, 1.0)
            || numelyx_matrix_create_dense(1, 1, huge_value, &huge)
                != NUMELYX_SUCCESS
            || numelyx_matrix_condition_estimate_one(huge, &scalar)
                != NUMELYX_SUCCESS || !close_to(scalar, 1.0)) {
            numelyx_matrix_release(huge);
            numelyx_matrix_release(tiny);
            goto fail;
        }

        numelyx_matrix_release(huge);
        numelyx_matrix_release(tiny);
    }
    {
        int positive_definite = 0, positive_semidefinite = 0;
        if (numelyx_matrix_classify_definiteness(
                a, &positive_definite, &positive_semidefinite)
                != NUMELYX_SUCCESS
            || !positive_definite || !positive_semidefinite)
            goto fail;
    }
    if (numelyx_matrix_sum(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 5.0)
        || numelyx_matrix_mean(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 1.25)
        || numelyx_matrix_variance(a, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 1.6875)
        || numelyx_matrix_variance(a, 1, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 2.25)
        || numelyx_matrix_variance(a, 2, &scalar)
            != NUMELYX_INVALID_ARGUMENT
        || numelyx_matrix_standard_deviation(a, 0, &scalar)
            != NUMELYX_SUCCESS || !close_to(scalar, sqrt(1.6875))
        || numelyx_matrix_standard_deviation(a, 1, &scalar)
            != NUMELYX_SUCCESS || !close_to(scalar, 1.5)
        || numelyx_matrix_standard_deviation(a, 2, &scalar)
            != NUMELYX_INVALID_ARGUMENT
        || numelyx_matrix_min(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 0.0)
        || numelyx_matrix_max(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 3.0))
        goto fail;

    numelyx_matrix_release(a); a = NULL;
    numelyx_matrix_release(b); b = NULL;
    if (numelyx_matrix_create_dense(3, 2, design, &a) != NUMELYX_SUCCESS
        || numelyx_matrix_create_dense(3, 1, response, &b) != NUMELYX_SUCCESS
        || numelyx_matrix_least_squares(a, b, &solution) != NUMELYX_SUCCESS
        || numelyx_matrix_get(solution, 0, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 1.0)
        || numelyx_matrix_get(solution, 1, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 2.0)
        || numelyx_matrix_classify_system(
            a, b, &classification) != NUMELYX_SUCCESS
        || classification != NUMELYX_MATRIX_SOLUTION_UNIQUE)
        goto fail;
    {
        const double deficient_values[] = {1.0, 2.0, 2.0, 4.0};
        const double consistent_values[] = {3.0, 6.0};
        const double inconsistent_values[] = {3.0, 7.0};
        const double consistent_multi_values[] = {3.0, 6.0, 6.0, 12.0};
        const double inconsistent_multi_values[] = {3.0, 6.0, 6.0, 13.0};
        numelyx_matrix *deficient = NULL;
        numelyx_matrix *consistent = NULL;
        numelyx_matrix *inconsistent = NULL;
        numelyx_matrix *consistent_multi = NULL;
        numelyx_matrix *inconsistent_multi = NULL;

        if (numelyx_matrix_create_dense(2, 2, deficient_values, &deficient)
                != NUMELYX_SUCCESS
            || numelyx_matrix_create_dense(2, 1, consistent_values, &consistent)
                != NUMELYX_SUCCESS
            || numelyx_matrix_create_dense(2, 1, inconsistent_values, &inconsistent)
                != NUMELYX_SUCCESS
            || numelyx_matrix_create_dense(
                2, 2, consistent_multi_values, &consistent_multi)
                != NUMELYX_SUCCESS
            || numelyx_matrix_create_dense(
                2, 2, inconsistent_multi_values, &inconsistent_multi)
                != NUMELYX_SUCCESS
            || numelyx_matrix_classify_system(
                deficient, consistent, &classification) != NUMELYX_SUCCESS
            || classification != NUMELYX_MATRIX_SOLUTION_INFINITE
            || numelyx_matrix_classify_system(
                deficient, inconsistent, &classification) != NUMELYX_SUCCESS
            || classification != NUMELYX_MATRIX_SOLUTION_INCONSISTENT
            || numelyx_matrix_classify_system(
                deficient, consistent_multi, &classification) != NUMELYX_SUCCESS
            || classification != NUMELYX_MATRIX_SOLUTION_INFINITE
            || numelyx_matrix_classify_system(
                deficient, inconsistent_multi, &classification) != NUMELYX_SUCCESS
            || classification != NUMELYX_MATRIX_SOLUTION_INCONSISTENT) {
            numelyx_matrix_release(inconsistent_multi);
            numelyx_matrix_release(consistent_multi);
            numelyx_matrix_release(inconsistent);
            numelyx_matrix_release(consistent);
            numelyx_matrix_release(deficient);
            goto fail;
        }

        numelyx_matrix_release(inconsistent_multi);
        numelyx_matrix_release(consistent_multi);
        numelyx_matrix_release(inconsistent);
        numelyx_matrix_release(consistent);
        numelyx_matrix_release(deficient);
    }

    numelyx_matrix_release(solution); solution = NULL;
    if (numelyx_matrix_weighted_least_squares(
            a, b, weights, 3, &solution) != NUMELYX_SUCCESS
        || numelyx_matrix_get(solution, 0, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 1.0)
        || numelyx_matrix_get(solution, 1, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 2.0))
        goto fail;

    if (numelyx_matrix_create_dense(2, 2, tri, &triangular) != NUMELYX_SUCCESS
        || numelyx_matrix_create_dense(2, 1, rhs_values, &rhs2)
            != NUMELYX_SUCCESS
        || numelyx_matrix_solve_triangular(
            triangular, rhs2, 1, 0, &tri_solution) != NUMELYX_SUCCESS
        || numelyx_matrix_get(tri_solution, 0, 0, &scalar)
            != NUMELYX_SUCCESS || !close_to(scalar, 1.0)
        || numelyx_matrix_get(tri_solution, 1, 0, &scalar)
            != NUMELYX_SUCCESS || !close_to(scalar, 1.0))
        goto fail;

    numelyx_matrix_release(extra);
    numelyx_matrix_release(basis);
    numelyx_matrix_release(a);
    numelyx_matrix_release(b);
    numelyx_matrix_release(q);
    numelyx_matrix_release(r);
    numelyx_matrix_release(u);
    numelyx_matrix_release(s);
    numelyx_matrix_release(vt);
    numelyx_matrix_release(pinv);
    numelyx_matrix_release(lower);
    numelyx_matrix_release(solution);
    numelyx_matrix_release(triangular);
    numelyx_matrix_release(rhs2);
    numelyx_matrix_release(tri_solution);
    return 0;

fail:
    numelyx_matrix_release(extra);
    numelyx_matrix_release(basis);
    numelyx_matrix_release(a);
    numelyx_matrix_release(b);
    numelyx_matrix_release(q);
    numelyx_matrix_release(r);
    numelyx_matrix_release(u);
    numelyx_matrix_release(s);
    numelyx_matrix_release(vt);
    numelyx_matrix_release(pinv);
    numelyx_matrix_release(lower);
    numelyx_matrix_release(solution);
    numelyx_matrix_release(triangular);
    numelyx_matrix_release(rhs2);
    numelyx_matrix_release(tri_solution);
    return 1;
}
