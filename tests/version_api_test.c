#include <string.h>

#include <numelyx/version.h>
#include <numelyx/status.h>

int main(void)
{
    if (numelyx_core_version_major() != UINT32_C(0)) {
        return 1;
    }
    if (numelyx_core_version_minor() != UINT32_C(1)) {
        return 2;
    }
    if (numelyx_core_version_patch() != UINT32_C(0)) {
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
    if (strcmp(numelyx_status_string((numelyx_status) 999), "unknown status") != 0) {
        return 7;
    }

    return 0;
}
