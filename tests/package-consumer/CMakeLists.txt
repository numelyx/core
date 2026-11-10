cmake_minimum_required(VERSION 3.16)
project(NumelyxCorePackageConsumer LANGUAGES C)

find_package(NumelyxCore CONFIG REQUIRED)

add_executable(numelyx_package_smoke main.c)
target_link_libraries(numelyx_package_smoke PRIVATE Numelyx::Core)
