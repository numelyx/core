#ifndef NUMELYX_VERSION_H
#define NUMELYX_VERSION_H

#include <stdint.h>

#include <numelyx/export.h>

NUMELYX_BEGIN_DECLS

/**
 * Return the Core semantic version components.
 *
 * These values describe the built library, not the stability of its C ABI.
 */
NUMELYX_API uint32_t numelyx_core_version_major(void);
NUMELYX_API uint32_t numelyx_core_version_minor(void);
NUMELYX_API uint32_t numelyx_core_version_patch(void);

/**
 * Return the supported ABI version.
 *
 * Zero means that the public C ABI is not yet declared stable. A non-zero
 * value will identify a compatibility contract and may change only through
 * the documented ABI-versioning policy.
 */
NUMELYX_API uint32_t numelyx_core_abi_version(void);

NUMELYX_END_DECLS

#endif /* NUMELYX_VERSION_H */
