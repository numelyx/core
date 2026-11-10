# Numelyx Core

Numelyx Core is a standalone C numerical-computing library migrated from the
native computational layer of [Numerus](https://github.com/ephpicman/Numerus).

The repository contains the language-independent numerical implementation:
immutable matrix storage, composable Matrix primitives and views, numerical
linear algebra, random-number generation, probability distributions, and
optimization primitives. PHP extension glue, Zend-specific public APIs,
PHPT tests, and PHP build configuration are intentionally excluded.

## Build

Requirements: CMake 3.16+ and a C11 compiler.

Static library (default):

```sh
cmake -S . -B build -DNUMELYX_BUILD_TESTS=ON
cmake --build build --config Debug --parallel
ctest --test-dir build -C Debug --output-on-failure
```

Shared library (for runtime loading and language-binding integration):

```sh
cmake -S . -B build-shared -DNUMELYX_BUILD_SHARED=ON -DNUMELYX_BUILD_TESTS=ON
cmake --build build-shared --config Debug --parallel
ctest --test-dir build-shared -C Debug --output-on-failure
```

The shared build produces the platform-native library (`.so`, `.dylib`, or
`.dll`). Windows builds currently use CMake's automatic symbol export. This
does not yet constitute a stable public C ABI: a curated public-header layer,
explicit API/export macros, and documented ownership/error contracts are
still required before external SDKs should bind to it. Internal headers and
symbol names remain subject to change.

Standalone builds use the C runtime allocator. The source retains conditional
Zend allocator support only in the upstream integration points; the CMake
target selects the libc allocator paths and does not include PHP headers.

## Install and consume from another CMake project

Install the selected static or shared build:

```sh
cmake --install build --prefix install
```

Then consume the exported target:

```cmake
find_package(NumelyxCore CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE Numelyx::Core)
```

The first supported header is `<numelyx/version.h>`, which exposes semantic-version and ABI-version introspection through `NUMELYX_API` declarations. `numelyx_core_abi_version()` currently returns `0`, explicitly indicating that no stable ABI has been declared.

The current transitional umbrella header is `<numelyx/numelyx.h>`. It includes
the existing subsystem headers installed under `numelyx/detail/`. These are
currently **unstable integration headers**, not a frozen public API: names,
layouts, and contracts may change until the dedicated public C API/ABI review
is complete. SDKs should not yet promise compatibility across Core releases.

## Sanitizer builds

The sanitizer configurations are intended for Linux with GCC or Clang and
cannot be enabled together.

AddressSanitizer, LeakSanitizer, and UndefinedBehaviorSanitizer:

```sh
cmake -S . -B build-asan -DNUMELYX_BUILD_TESTS=ON -DNUMELYX_ENABLE_ASAN_UBSAN=ON
cmake --build build-asan --config Debug --parallel
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 \
    ctest --test-dir build-asan -C Debug --output-on-failure
```

ThreadSanitizer:

```sh
cmake -S . -B build-tsan -DNUMELYX_BUILD_TESTS=ON -DNUMELYX_ENABLE_TSAN=ON
cmake --build build-tsan --config Debug --parallel
TSAN_OPTIONS=halt_on_error=1 ctest --test-dir build-tsan -C Debug --output-on-failure
```

The allocation-failure tests that use GNU linker `--wrap` are built only with
GNU on Linux. Other compiler/platform configurations skip those tests.

## Thread safety

Shared opaque objects use atomic reference counting. Matrix lazy caches are
serialized per Matrix, and mutable RNG state is serialized per RNG. See
[`docs/thread-safety-and-lifetime.md`](docs/thread-safety-and-lifetime.md) for
ownership requirements and callback restrictions. The planned public boundary is documented in [`docs/public-api-and-abi.md`](docs/public-api-and-abi.md).

## Benchmarks

```sh
cmake -S . -B build -DNUMELYX_BUILD_BENCHMARKS=ON
cmake --build build --target matrix_benchmark
./build/matrix_benchmark
```

## Source layout

- `src/core`: checked size arithmetic, stable scalar functions, generic object lifetime, and synchronization primitives.
- `src/storage`: opaque immutable storage representations.
- `src/matrix`: composable Matrix primitives, views, and factorizations.
- `src/statistics`: RNG, probability, and optimization primitives.
- `tests`: native C regression tests.
- `docs`: numerical, subsystem, thread-safety, lifetime, and public ABI contracts.

The C symbols retain their upstream `numelyx_core_*` names during the extraction;
renaming the ABI is a separate compatibility decision.
