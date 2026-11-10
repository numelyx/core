cmake_minimum_required(VERSION 3.16)

project(NumelyxCore VERSION 0.1.0 LANGUAGES C)

option(NUMELYX_BUILD_TESTS "Build native C tests" ON)
option(NUMELYX_BUILD_BENCHMARKS "Build diagnostic benchmarks" OFF)
option(NUMELYX_BUILD_SHARED "Build Numelyx Core as a shared library" OFF)
option(NUMELYX_ENABLE_TSAN "Enable ThreadSanitizer (Linux, GCC/Clang)" OFF)
option(NUMELYX_ENABLE_ASAN_UBSAN "Enable AddressSanitizer and UndefinedBehaviorSanitizer (Linux, GCC/Clang)" OFF)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

file(GLOB_RECURSE NUMELYX_CORE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/storage/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/matrix/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/statistics/*.c"
)

if(NUMELYX_BUILD_SHARED)
    set(NUMELYX_LIBRARY_KIND SHARED)
else()
    set(NUMELYX_LIBRARY_KIND STATIC)
endif()

add_library(numelyx_core ${NUMELYX_LIBRARY_KIND} ${NUMELYX_CORE_SOURCES})
add_library(Numelyx::Core ALIAS numelyx_core)

# Shared builds must be usable by downstream language bindings. The explicit
# export macro covers the supported version API; detail symbols remain visible
# during the transition to the curated public API.
set_target_properties(numelyx_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

find_package(Threads REQUIRED)
target_link_libraries(numelyx_core PUBLIC Threads::Threads)
target_compile_definitions(numelyx_core
    PRIVATE
        NUMELYX_BUILDING_LIBRARY
        NUMELYX_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
        NUMELYX_VERSION_MINOR=${PROJECT_VERSION_MINOR}
        NUMELYX_VERSION_PATCH=${PROJECT_VERSION_PATCH}
)
if(NUMELYX_BUILD_SHARED)
    target_compile_definitions(numelyx_core PUBLIC NUMELYX_SHARED)
endif()

if(NUMELYX_ENABLE_TSAN)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
        message(FATAL_ERROR "NUMELYX_ENABLE_TSAN is supported only on Linux")
    endif()
    if(NOT CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "NUMELYX_ENABLE_TSAN requires GCC or Clang")
    endif()

    target_compile_options(numelyx_core PUBLIC
        -fsanitize=thread
        -fno-omit-frame-pointer
        -g
    )
    target_link_options(numelyx_core PUBLIC -fsanitize=thread)
endif()

if(NUMELYX_ENABLE_TSAN AND NUMELYX_ENABLE_ASAN_UBSAN)
    message(FATAL_ERROR "ThreadSanitizer cannot be combined with Address/UndefinedBehaviorSanitizer")
endif()

if(NUMELYX_ENABLE_ASAN_UBSAN)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
        message(FATAL_ERROR "NUMELYX_ENABLE_ASAN_UBSAN is supported only on Linux")
    endif()
    if(NOT CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "NUMELYX_ENABLE_ASAN_UBSAN requires GCC or Clang")
    endif()

    target_compile_options(numelyx_core PUBLIC
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
        -fno-sanitize-recover=all
        -g
    )
    target_link_options(numelyx_core PUBLIC
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
    )
endif()

if(MSVC)
    target_compile_options(numelyx_core PUBLIC /experimental:c11atomics)
endif()

include(GNUInstallDirs)

target_include_directories(numelyx_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/storage>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/matrix>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/statistics>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/numelyx/detail/core>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/numelyx/detail/storage>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/numelyx/detail/matrix>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/numelyx/detail/statistics>
)

set_target_properties(numelyx_core PROPERTIES EXPORT_NAME Core)

# The current declarations are installed as explicitly unstable detail headers.
# They are useful for in-repository and early SDK integration, but are not yet
# the compatibility promise of a versioned public C API.
set(NUMELYX_GENERATED_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/include")
file(MAKE_DIRECTORY "${NUMELYX_GENERATED_INCLUDE_DIR}/numelyx")
set(NUMELYX_UMBRELLA_HEADER
    "/* Supported public API umbrella header. */\\n#ifndef NUMELYX_NUMELYX_H\\n#define NUMELYX_NUMELYX_H\\n#include <numelyx/version.h>\\n#include <numelyx/status.h>\\n#include <numelyx/matrix.h>\\n#include <numelyx/rng.h>\\n#include <numelyx/optimizer.h>\\n#include <numelyx/probability.h>\\n#endif /* NUMELYX_NUMELYX_H */\\n"
)
file(WRITE "${NUMELYX_GENERATED_INCLUDE_DIR}/numelyx/numelyx.h" "${NUMELYX_UMBRELLA_HEADER}")

target_include_directories(numelyx_core
    PUBLIC
        $<BUILD_INTERFACE:${NUMELYX_GENERATED_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS numelyx_core
    EXPORT NumelyxCoreTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

foreach(subsystem IN ITEMS core storage matrix statistics)
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/${subsystem}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/numelyx/detail/${subsystem}"
        FILES_MATCHING PATTERN "*.h"
    )
endforeach()

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(FILES "${NUMELYX_GENERATED_INCLUDE_DIR}/numelyx/numelyx.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/numelyx"
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/NumelyxCoreConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/NumelyxCoreConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/NumelyxCore"
)
if(PROJECT_VERSION_MAJOR EQUAL 0)
    # In pre-1.0 releases, minor-version changes may be breaking.
    set(NUMELYX_PACKAGE_COMPATIBILITY SameMinorVersion)
else()
    set(NUMELYX_PACKAGE_COMPATIBILITY SameMajorVersion)
endif()

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/NumelyxCoreConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY ${NUMELYX_PACKAGE_COMPATIBILITY}
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/NumelyxCoreConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/NumelyxCoreConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/NumelyxCore"
)
install(EXPORT NumelyxCoreTargets
    FILE NumelyxCoreTargets.cmake
    NAMESPACE Numelyx::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/NumelyxCore"
)

if(UNIX)
    target_link_libraries(numelyx_core PUBLIC m)
endif()

if(NUMELYX_BUILD_TESTS)
    enable_testing()
    file(GLOB NUMELYX_TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*_test.c"
    )

    set(NUMELYX_ALLOC_WRAP_TESTS
        matrix_allocation_failure_test
        matrix_concurrent_allocation_failure_test
        storage_allocation_failure_test
    )

    foreach(test_source IN LISTS NUMELYX_TEST_SOURCES)
        get_filename_component(test_name "${test_source}" NAME_WE)

        # Matrix and Storage allocation injection uses GNU ld --wrap.
        if(test_name IN_LIST NUMELYX_ALLOC_WRAP_TESTS AND NUMELYX_BUILD_SHARED)
            message(STATUS "Skipping ${test_name}: GNU ld --wrap cannot inject allocations into a shared library")
            continue()
        endif()

        if(test_name IN_LIST NUMELYX_ALLOC_WRAP_TESTS
           AND NOT (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND UNIX AND NOT APPLE))
            message(STATUS "Skipping ${test_name}: GNU ld --wrap is required")
            continue()
        endif()

        add_executable(${test_name} "${test_source}")
        target_link_libraries(${test_name} PRIVATE Numelyx::Core)

        if(test_name IN_LIST NUMELYX_ALLOC_WRAP_TESTS)
            target_link_options(${test_name} PRIVATE
                "-Wl,--wrap=malloc"
                "-Wl,--wrap=calloc"
                "-Wl,--wrap=free"
            )
        elseif(test_name MATCHES "^(gaussian_log_density_test|gaussian_sampling_test)$")
            target_sources(${test_name} PRIVATE src/statistics/numelyx_core_probability.c)
            target_compile_definitions(${test_name} PRIVATE NUMERUS_PROBABILITY_TEST_ALLOCATOR)
        elseif(test_name STREQUAL "optimizer_test")
            target_sources(${test_name} PRIVATE src/statistics/numelyx_core_optimizer.c)
            target_compile_definitions(${test_name} PRIVATE NUMERUS_OPTIMIZER_TEST_ALLOCATOR)
        elseif(test_name MATCHES "^(rng_sampling_test|rng_state_test)$")
            target_sources(${test_name} PRIVATE
                src/statistics/numelyx_core_rng.c
                src/statistics/numelyx_core_rng_sampling.c
            )
            target_compile_definitions(${test_name} PRIVATE NUMERUS_RNG_TEST_ALLOCATOR)
        endif()

        if(test_name STREQUAL "lifetime_test")
            target_link_libraries(${test_name} PRIVATE Threads::Threads)
        endif()

        add_test(NAME ${test_name} COMMAND ${test_name})
        set_tests_properties(${test_name} PROPERTIES TIMEOUT 60)
    endforeach()
endif()

if(NUMELYX_BUILD_BENCHMARKS)
    add_executable(matrix_benchmark benchmarks/matrix_benchmark.c)
    target_link_libraries(matrix_benchmark PRIVATE Numelyx::Core)
    if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND UNIX AND NOT APPLE)
        target_link_options(matrix_benchmark PRIVATE
            "-Wl,--wrap=malloc"
            "-Wl,--wrap=calloc"
            "-Wl,--wrap=free"
        )
    endif()
endif()
