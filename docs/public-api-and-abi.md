# Public C API and ABI Policy

## Status

This document defines the boundary required before Numelyx Core is consumed as
a supported dependency by independently released language SDKs. It does not
declare the current `src/` headers or exported symbols stable.

The currently installed headers under `numelyx/detail/` are transitional
integration headers. They expose implementation-era declarations and may
change without ABI compatibility guarantees.

## Public boundary

The supported API will live under `include/numelyx/`. Consumers must not need
headers from `src/`, implementation structs, allocator macros, mutex types,
or subsystem-private status definitions.

Public objects are opaque handles. Storage layout, cache state, reference
counters, synchronization primitives, and factorization internals remain
private. Public structs are limited to versioned value/configuration types
whose layout is intentionally part of the ABI.

## ABI requirements

Before declaring the C ABI stable, the implementation must provide:

1. A platform-neutral `NUMELYX_API` import/export macro and C++ linkage guards.
2. A documented ABI version query independent of the package's semantic version.
3. A single documented public status-code namespace, with subsystem-specific
   details represented by stable error values or documented extended status.
4. Explicit ownership rules for every returned handle and allocated buffer.
   Every allocation crossing the ABI must have a matching Core release function;
   callers must not free Core memory with an unrelated allocator.
5. Opaque handles and fixed-width types where the ABI requires a fixed width.
   `size_t` is reserved for process-local dimensions and buffer lengths.
6. Defined behavior for NULL arguments, output initialization on failure,
   overflow, non-finite inputs, and callback failures.
7. Documented thread-safety and callback-context lifetime rules.
8. Symbol visibility tests and C consumer tests on Linux, Windows, and macOS.

## Compatibility policy

- Internal `src/` headers and symbols are not API.
- The transitional umbrella header is not a compatibility promise.
- Public declarations are introduced deliberately; internal names are not
  mechanically renamed until their public replacements and migration path exist.
- ABI-breaking changes require an ABI version change and an explicit release
  decision.
- Static and shared builds must expose the same source-level public API.
- The Core library must not depend on PHP, Zend Engine, HTTP, JSON, or a server
  transport. Those belong in adapters and SDKs.

## Initial API surface

The first supported surface should be intentionally small:

- Core version and ABI information.
- Status/error inspection.
- Matrix creation from supported storage inputs, retain/release, dimensions,
  element access, and views.
- Explicitly named numerical operations added in later API increments.
- RNG creation, deterministic state operations, sampling, and release.

Do not export every internal helper simply because it exists. The API should
expose composable numerical primitives rather than model-specific workflows.
