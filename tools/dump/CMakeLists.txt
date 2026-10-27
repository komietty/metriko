# metriko_dump: runs compute_quadrangulation_impl step by step and writes every stage as .npy (see FORMAT.md).
# built the way the reference is: Eigen's solvers (no SuiteSparse), no OpenMP, serial libigl, no FMA contraction.
#
#   cmake -S tools/dump -B tools/dump/build \
#         -DFETCHCONTENT_SOURCE_DIR_LIBIGL=<libigl v2.6.4> -DFETCHCONTENT_SOURCE_DIR_EIGEN=<eigen 5.0.1>
#   cmake --build tools/dump/build -j8
cmake_minimum_required(VERSION 3.28)
project(metriko_dump CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG" CACHE STRING "" FORCE)

set(METRIKO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH "Path to the metriko repository")
set(SUITESPARSE   OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
add_subdirectory(${METRIKO_DIR} metriko)

add_executable(metriko_dump dump.cpp)
target_link_libraries(metriko_dump PRIVATE metriko)
target_compile_definitions(metriko_dump PRIVATE IGL_PARALLEL_FOR_FORCE_SERIAL)
target_compile_options(metriko_dump PRIVATE -ffp-contract=off)
