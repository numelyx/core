#include <numelyx/numelyx.h>

int main(void)
{
    numelyx_core_rng *rng = NULL;

    if (numelyx_core_rng_create(UINT64_C(42), UINT64_C(54), &rng)
        != NUMERUS_RNG_SUCCESS) {
        return 1;
    }

    numelyx_core_rng_destroy(rng);
    return 0;
}
