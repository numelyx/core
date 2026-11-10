#ifndef NUMELYX_STATUS_H
#define NUMELYX_STATUS_H

#include <numelyx/export.h>

NUMELYX_BEGIN_DECLS

/** Public names for errors exposed by the supported C API. */
typedef enum {
    NUMELYX_SUCCESS = 0,
    NUMELYX_INVALID_ARGUMENT,
    NUMELYX_OUT_OF_MEMORY,
    NUMELYX_SIZE_OVERFLOW,
    NUMELYX_OUT_OF_BOUNDS,
    NUMELYX_DIMENSION_MISMATCH,
    NUMELYX_NON_FINITE,
    NUMELYX_SINGULAR,
    NUMELYX_RANK_DEFICIENT,
    NUMELYX_NO_CONVERGENCE,
    NUMELYX_DIVISION_BY_ZERO,
    NUMELYX_NOT_SYMMETRIC,
    NUMELYX_NOT_POSITIVE_DEFINITE,
    NUMELYX_NOT_POSITIVE_SEMIDEFINITE,
    NUMELYX_PIVOT_TOO_SMALL,
    NUMELYX_NO_OBJECTIVE_VALUE,
    NUMELYX_NUMERICAL_ERROR
} numelyx_status;

/** Return a static, NUL-terminated description for a public status value. */
NUMELYX_API const char *numelyx_status_string(numelyx_status status);

NUMELYX_END_DECLS

#endif /* NUMELYX_STATUS_H */
