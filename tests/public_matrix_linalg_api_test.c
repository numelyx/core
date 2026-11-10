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
    size_t rank = 0;
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

    if (numelyx_matrix_symmetric_eigen(
            a, &u, &vt) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(u) != 2
        || numelyx_matrix_columns(u) != 1
        || numelyx_matrix_rows(vt) != 2
        || numelyx_matrix_columns(vt) != 2)
        goto fail;
    numelyx_matrix_release(u); u = NULL;
    numelyx_matrix_release(vt); vt = NULL;
    if (numelyx_matrix_ldlt_decompose(
            a, &u, &vt) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(u) != 2
        || numelyx_matrix_rows(vt) != 2)
        goto fail;
    numelyx_matrix_release(u); u = NULL;
    numelyx_matrix_release(vt); vt = NULL;
    if (numelyx_matrix_spectral_norm(a, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 2.0)
        || numelyx_matrix_condition_estimate_one(a, &scalar)
            != NUMELYX_SUCCESS || !close_to(scalar, 2.0))
        goto fail;
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
        || !close_to(scalar, 1.1875)
        || numelyx_matrix_variance(a, 2, &scalar)
            != NUMELYX_INVALID_ARGUMENT)
        goto fail;

    numelyx_matrix_release(a); a = NULL;
    numelyx_matrix_release(b); b = NULL;
    if (numelyx_matrix_create_dense(3, 2, design, &a) != NUMELYX_SUCCESS
        || numelyx_matrix_create_dense(3, 1, response, &b) != NUMELYX_SUCCESS
        || numelyx_matrix_least_squares(a, b, &solution) != NUMELYX_SUCCESS
        || numelyx_matrix_get(solution, 0, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 1.0)
        || numelyx_matrix_get(solution, 1, 0, &scalar) != NUMELYX_SUCCESS
        || !close_to(scalar, 2.0))
        goto fail;

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
