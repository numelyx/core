#include <numelyx/version.h>

#ifndef NUMELYX_VERSION_MAJOR
#    define NUMELYX_VERSION_MAJOR 0
#endif
#ifndef NUMELYX_VERSION_MINOR
#    define NUMELYX_VERSION_MINOR 1
#endif
#ifndef NUMELYX_VERSION_PATCH
#    define NUMELYX_VERSION_PATCH 0
#endif

uint32_t numelyx_core_version_major(void)
{
    return (uint32_t) NUMELYX_VERSION_MAJOR;
}

uint32_t numelyx_core_version_minor(void)
{
    return (uint32_t) NUMELYX_VERSION_MINOR;
}

uint32_t numelyx_core_version_patch(void)
{
    return (uint32_t) NUMELYX_VERSION_PATCH;
}

uint32_t numelyx_core_abi_version(void)
{
    /* The public API is still transitional; no stable ABI has been declared. */
    return UINT32_C(0);
}
