#include <numelyx/numelyx.h>

int main(void)
{
    numelyx_core_rng *rng = NULL;

    if (numelyx_core_version_major() != UINT32_C(0)
        || numelyx_core_version_minor() != UINT32_C(1)
        || numelyx_core_version_patch() != UINT32_C(0)
        || numelyx_core_abi_version() != UINT32_C(0)) {
        return 1;
    }

    if (numelyx_core_rng_create(UINT64_C(42), UINT64_C(54), &rng)
        != NUMERUS_RNG_SUCCESS) {
        return 2;
    }

    numelyx_core_rng_destroy(rng);
    return 0;
}
