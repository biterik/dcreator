cmake_minimum_required(VERSION 3.20)
project(dcreator
    VERSION 2.0.0
    DESCRIPTION "Parallel dislocation inserter for atomistic simulations"
    LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(DCREATOR_ENABLE_OPENMP "Enable OpenMP parallelization" ON)
option(DCREATOR_BUILD_TESTS   "Build regression and unit tests" ON)
option(DCREATOR_NATIVE_ARCH   "Compile with -march=native" OFF)

# -------------------------------------------------------------------
# macOS: prefer Homebrew LLVM for OpenMP support. Apple clang ships
# without libomp and needs fragile workarounds; LLVM clang just works.
# -------------------------------------------------------------------
if(APPLE AND NOT DEFINED ENV{DCREATOR_NO_HOMEBREW_LLVM})
    foreach(_prefix "/opt/homebrew/opt/llvm" "/usr/local/opt/llvm")
        if(EXISTS "${_prefix}/bin/clang")
            message(STATUS "Using Homebrew LLVM at ${_prefix}")
            list(APPEND CMAKE_PREFIX_PATH "${_prefix}")
            if(NOT CMAKE_C_COMPILER)
                set(CMAKE_C_COMPILER "${_prefix}/bin/clang" CACHE FILEPATH "" FORCE)
            endif()
            break()
        endif()
    endforeach()
    foreach(_libomp "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")
        if(EXISTS "${_libomp}")
            list(APPEND CMAKE_PREFIX_PATH "${_libomp}")
            break()
        endif()
    endforeach()
endif()

if(DCREATOR_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS C)
endif()

# -------------------------------------------------------------------
# libdcreator
# -------------------------------------------------------------------
add_library(dcreator_core STATIC
    src/core/geometry.c
    src/core/kernel.c
    src/core/params.c
    src/core/atoms.c
    src/io/format_detect.c
    src/io/lammps_data.c
    src/io/lammps_dump.c
    src/io/imd.c
    src/io/io_util.c
)
target_include_directories(dcreator_core PUBLIC include)
target_compile_options(dcreator_core PRIVATE
    -Wall -Wextra -Wpedantic -Wshadow -Wformat=2
    -Wno-unused-parameter
    -Wno-gnu-zero-variadic-macro-arguments
)
if(DCREATOR_NATIVE_ARCH)
    target_compile_options(dcreator_core PRIVATE -march=native)
endif()
target_link_libraries(dcreator_core PUBLIC m)
if(DCREATOR_ENABLE_OPENMP)
    target_link_libraries(dcreator_core PUBLIC OpenMP::OpenMP_C)
    target_compile_definitions(dcreator_core PUBLIC DCREATOR_HAVE_OPENMP=1)
endif()

# -------------------------------------------------------------------
# dcreator CLI
# -------------------------------------------------------------------
add_executable(dcreator src/cli/main.c)
target_link_libraries(dcreator PRIVATE dcreator_core)

install(TARGETS dcreator RUNTIME DESTINATION bin)
install(TARGETS dcreator_core ARCHIVE DESTINATION lib)
install(FILES include/dcreator.h DESTINATION include)

# -------------------------------------------------------------------
# Tests
# -------------------------------------------------------------------
if(DCREATOR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
