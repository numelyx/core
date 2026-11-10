#include <math.h>
#include <stddef.h>

#include <numelyx/matrix.h>

static int close_to(double actual, double expected)
{
    return fabs(actual - expected) < 1e-8;
}

int main(void)
{
    numelyx_matrix *a = NULL;
    numelyx_matrix *b = NULL;
    numelyx_matrix *q = NULL;
    numelyx_matrix *r = NULL;
    numelyx_matrix *u = NULL;
    numelyx_matrix *s = NULL;
    numelyx_matrix *vt = NULL;
    numelyx_matrix *pinv = NULL;
    numelyx_matrix *lower = NULL;
    numelyx_matrix *solution = NULL;
    numelyx_matrix *triangular = NULL;
    size_t rank = 0;
    int sign = 0;
    double scalar = 0.0;
    const double diagonal[] = {2.0, 0.0, 0.0, 3.0};
    const double design[] = {1.0, 0.0, 1.0, 1.0, 1.0, 2.0};
    const double response[] = {1.0, 3.0, 5.0};
    const double weights[] = {1.0, 2.0, 1.0};
    const double spd[] = {4.0, 2.0, 2.0, 3.0};
    const double tri[] = {2.0, 1.0, 0.0, 3.0};
    double log_abs = 0.0;

    if (numelyx_matrix_create_dense(2, 2, diagonal, &a) != NUMELYX_SUCCESS)
        return 1;
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

    if (numelyx_matrix_qr_decompose(a, &q, &r) != NUMELYX_SUCCESS
        || numelyx_matrix_rows(q) != 2 || numelyx_matrix_columns(q) != 2
        || numelyx_matrix_rows(r) != 2 || numelyx_matrix_columns(r) != 2)
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

    numelyx_matrix_release(a); a = NULL;
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
        || numelyx_matrix_solve_triangular(
            triangular, b, 1, 0, &q) == NUMELYX_SUCCESS)
        goto fail; /* RHS has the wrong row count. */

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
    return 2;
}
