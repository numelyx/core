#include <numelyx/status.h>

const char *numelyx_status_string(numelyx_status status)
{
    switch (status) {
        case NUMELYX_SUCCESS:
            return "success";
        case NUMELYX_INVALID_ARGUMENT:
            return "invalid argument";
        case NUMELYX_OUT_OF_MEMORY:
            return "out of memory";
        case NUMELYX_SIZE_OVERFLOW:
            return "size overflow";
        case NUMELYX_OUT_OF_BOUNDS:
            return "out of bounds";
        case NUMELYX_DIMENSION_MISMATCH:
            return "dimension mismatch";
        case NUMELYX_NON_FINITE:
            return "non-finite value";
        case NUMELYX_SINGULAR:
            return "singular matrix";
        case NUMELYX_RANK_DEFICIENT:
            return "rank deficient";
        case NUMELYX_NO_CONVERGENCE:
            return "no convergence";
        case NUMELYX_DIVISION_BY_ZERO:
            return "division by zero";
        case NUMELYX_NOT_SYMMETRIC:
            return "not symmetric";
        case NUMELYX_NOT_POSITIVE_DEFINITE:
            return "not positive definite";
        case NUMELYX_NOT_POSITIVE_SEMIDEFINITE:
            return "not positive semidefinite";
        case NUMELYX_PIVOT_TOO_SMALL:
            return "pivot too small";
        case NUMELYX_NO_OBJECTIVE_VALUE:
            return "no objective value";
        case NUMELYX_NUMERICAL_ERROR:
            return "numerical error";
        default:
            return "unknown status";
    }
}
