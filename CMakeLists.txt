cmake_minimum_required(VERSION 3.28)
project(metriko)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_OPENMP "Enable OpenMP support" ON)

#========= libigl (also provides Eigen) ================================================================================
if (NOT TARGET igl::core)
    include(FetchContent)
    FetchContent_Declare(
            libigl
            GIT_REPOSITORY https://github.com/libigl/libigl.git
            GIT_TAG v2.6.4
    )
    FetchContent_MakeAvailable(libigl)
endif()

# Work with non-standard homebrew installations
if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    find_program(HOMEBREW_EXECUTABLE brew)
    mark_as_advanced(FORCE HOMEBREW_EXECUTABLE)
    if (HOMEBREW_EXECUTABLE)
        # Detected a Homebrew install, query for its install prefix.
        execute_process(COMMAND ${HOMEBREW_EXECUTABLE} --prefix
                OUTPUT_VARIABLE HOMEBREW_INSTALL_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Detected Homebrew with install prefix: "
                "${HOMEBREW_INSTALL_PREFIX}, adding to CMake search paths.")
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_INSTALL_PREFIX}")
    endif()
endif()

### Handle windows-specific fixes
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    add_definitions (-DNOMINMAX)  # don't use weird windows built-in min/max
    add_definitions (-D_USE_MATH_DEFINES)  # match unix behavior of constants in cmath
endif()

#========= SuiteSparse (CHOLMOD / SPQR, via the CMake config shipped with SuiteSparse >= 7) ============================
set(SUITESPARSE_HELP "Enable SuiteSparse (CHOLMOD / SPQR). Falls back to Eigen's solvers when off or not found.")
option(SUITESPARSE "${SUITESPARSE_HELP}" ON)
set(GC_HAVE_SUITESPARSE FALSE)
if (SUITESPARSE)
    find_package(CHOLMOD CONFIG QUIET)
    find_package(SPQR    CONFIG QUIET)
    if (CHOLMOD_FOUND AND SPQR_FOUND)
        message(STATUS "Found SuiteSparse (CHOLMOD ${CHOLMOD_VERSION}, SPQR ${SPQR_VERSION}), building with SuiteSparse.")
        set(GC_HAVE_SUITESPARSE TRUE)
    else()
        message(STATUS "SuiteSparse (CHOLMOD / SPQR) not found, falling back to Eigen's solvers.")
        set(SUITESPARSE OFF CACHE BOOL "${SUITESPARSE_HELP}" FORCE)
    endif()
else()
    message(STATUS "Building without SuiteSparse, using Eigen's solvers.")
endif()

#========= Try finding OpenMp ==========================================================================================

set(METRIKO_HAVE_OPENMP FALSE)
if (ENABLE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        set(METRIKO_HAVE_OPENMP TRUE)
        message(STATUS "OpenMP found. CXX flags: ${OpenMP_CXX_FLAGS}")
    else()
        message(WARNING "OpenMP not found, building without OpenMP support.")
    endif()
endif()
#========== Done finding OpenMp ========================================================================================

add_library(${PROJECT_NAME} STATIC
        include/metriko/tmesh/emesh_collapse_range.cpp
        include/metriko/tmesh/emesh_collapse_thalf.cpp
        include/metriko/tmesh/emesh_collapse_tq_prepare.cpp
        include/metriko/tmesh/emesh_collapse_tq_execute.cpp
        include/metriko/tmesh/emesh_postprocess.cpp
)

#target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX _USE_MATH_DEFINES METRIKO_DEBUG)

target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)
target_include_directories(${PROJECT_NAME} PUBLIC include)

if(GC_HAVE_SUITESPARSE)
    target_compile_definitions(${PROJECT_NAME} INTERFACE GC_HAVE_SUITESPARSE)
    target_link_libraries(${PROJECT_NAME} INTERFACE SuiteSparse::CHOLMOD SuiteSparse::SPQR)
endif()

if(METRIKO_HAVE_OPENMP)
    target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)
endif()
