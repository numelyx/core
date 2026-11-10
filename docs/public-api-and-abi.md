# Public C API and ABI Policy

## Status

This document defines the boundary required before Numelyx Core is consumed as
a supported dependency by independently released language SDKs. It does not
declare the current `src/` headers or exported symbols stable.

Internal headers remain under `src/` and are not installed or exported to
consumers. They expose implementation-era declarations and may change without
ABI compatibility guarantees.

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
- The supported umbrella header groups only curated public headers, but it is
  not yet a compatibility promise while the ABI version remains zero.
- Public declarations are introduced deliberately; internal names are not
  mechanically renamed until their public replacements and migration path exist.
- ABI-breaking changes require an ABI version change and an explicit release
  decision.
- Static and shared builds must expose the same source-level public API.
- The Core library must not depend on PHP, Zend Engine, HTTP, JSON, or a server
  transport. Those belong in adapters and SDKs.

## Initial API surface

The first supported surface should be intentionally small:

- Core semantic version and ABI-stability introspection.
- A public status-code namespace with static status-description strings.
- Dense, identity, zero, constant, scaled-identity, diagonal, packed triangular, symmetric, default-value sparse, and fixed-bandwidth Matrix
  creation, retain/release,
  dimensions, element access, lazy transpose/slice/selection/flip/reshape/flatten,
  rotation, removal, swap, join, padding, repetition, block-grid, diagonal,
  and permutation views; lazy element-wise addition,
  subtraction, Hadamard product and scalar scaling, eager matrix multiplication, determinant, LU-based inverse and solve, and
  explicit materialization through opaque handles. Structured constructors preserve implicit storage rather than
  materializing dense data; views retain their parent dependencies.
- Public Matrix factorization and analysis primitives: QR, reduced SVD, Moore–Penrose pseudoinverse, Cholesky, LU solve, triangular solve, ordinary and weighted least squares, numerical rank, signed log-determinant, trace, and Frobenius/induced norms. These are thin status/ownership adapters over the Core algorithms; they do not introduce model-specific workflows.
- Public Matrix introspection exposes stable storage-kind identifiers and logical structural flags. A derived view reports the root storage kind while structural flags describe the view's logical values.
- RNG creation, exact-state cloning, retain/release, deterministic raw output,
  uniform/normal sampling, and index sampling with an explicit Core deallocator.
- BFGS optimizer options with explicit structure-size/version fields,
  objective/gradient callbacks, minimization, immutable result ownership, and
  result inspection through an opaque handle.
- Multivariate Gaussian log density and sampling using the existing Matrix and
  RNG handles; density uses Cholesky/log-determinant operations rather than
  explicit covariance inversion.

Do not export every internal helper simply because it exists. The API should
expose composable numerical primitives rather than model-specific workflows.
