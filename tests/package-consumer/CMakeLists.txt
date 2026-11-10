cmake_minimum_required(VERSION 3.16)
project(NumelyxCorePackageConsumer LANGUAGES C)

enable_testing()
find_package(NumelyxCore CONFIG REQUIRED)

add_executable(numelyx_package_smoke main.c)
target_link_libraries(numelyx_package_smoke PRIVATE Numelyx::Core)
target_compile_definitions(numelyx_package_smoke PRIVATE
    NUMELYX_EXPECTED_VERSION_MAJOR=${NumelyxCore_VERSION_MAJOR}
    NUMELYX_EXPECTED_VERSION_MINOR=${NumelyxCore_VERSION_MINOR}
    NUMELYX_EXPECTED_VERSION_PATCH=${NumelyxCore_VERSION_PATCH}
)
add_test(NAME numelyx_package_smoke COMMAND numelyx_package_smoke)
