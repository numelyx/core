#include <string.h>

#include <numelyx/version.h>
#include <numelyx/status.h>

#ifndef NUMELYX_EXPECTED_VERSION_MAJOR
#    define NUMELYX_EXPECTED_VERSION_MAJOR 0
#endif
#ifndef NUMELYX_EXPECTED_VERSION_MINOR
#    define NUMELYX_EXPECTED_VERSION_MINOR 1
#endif
#ifndef NUMELYX_EXPECTED_VERSION_PATCH
#    define NUMELYX_EXPECTED_VERSION_PATCH 0
#endif

int main(void)
{
    if (numelyx_core_version_major() != (uint32_t) NUMELYX_EXPECTED_VERSION_MAJOR) {
        return 1;
    }
    if (numelyx_core_version_minor() != (uint32_t) NUMELYX_EXPECTED_VERSION_MINOR) {
        return 2;
    }
    if (numelyx_core_version_patch() != (uint32_t) NUMELYX_EXPECTED_VERSION_PATCH) {
        return 3;
    }
    if (numelyx_core_abi_version() != UINT32_C(0)) {
        return 4;
    }
    if (strcmp(numelyx_status_string(NUMELYX_SUCCESS), "success") != 0) {
        return 5;
    }
    if (strcmp(numelyx_status_string(NUMELYX_NOT_POSITIVE_DEFINITE),
               "not positive definite") != 0) {
        return 6;
    }
    if (strcmp(numelyx_status_string(NUMELYX_NO_OBJECTIVE_VALUE),
               "no objective value") != 0) {
        return 7;
    }
    if (strcmp(numelyx_status_string((numelyx_status) 999), "unknown status") != 0) {
        return 8;
    }

    return 0;
}
