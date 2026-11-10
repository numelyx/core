#include <numelyx/version.h>

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

    return 0;
}
