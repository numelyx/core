cmake_minimum_required(VERSION 3.16)
project(NumelyxCorePackageConsumer LANGUAGES C)

enable_testing()
find_package(NumelyxCore CONFIG REQUIRED)

foreach(_property
    IMPORTED_CONFIGURATIONS
    IMPORTED_LOCATION
    IMPORTED_LOCATION_NOCONFIG
    IMPORTED_LOCATION_DEBUG
    IMPORTED_LOCATION_RELEASE
    IMPORTED_IMPLIB
    IMPORTED_IMPLIB_NOCONFIG
    IMPORTED_IMPLIB_DEBUG
    IMPORTED_IMPLIB_RELEASE
)
    get_target_property(_value Numelyx::Core ${_property})
    message(STATUS "Numelyx::Core ${_property}=${_value}")
endforeach()

add_executable(numelyx_package_smoke main.c)
target_link_libraries(numelyx_package_smoke PRIVATE Numelyx::Core)
target_compile_definitions(numelyx_package_smoke PRIVATE
    NUMELYX_EXPECTED_VERSION_MAJOR=${NumelyxCore_VERSION_MAJOR}
    NUMELYX_EXPECTED_VERSION_MINOR=${NumelyxCore_VERSION_MINOR}
    NUMELYX_EXPECTED_VERSION_PATCH=${NumelyxCore_VERSION_PATCH}
)
add_test(NAME numelyx_package_smoke COMMAND numelyx_package_smoke)
