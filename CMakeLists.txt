cmake_minimum_required(VERSION 3.28)
project(metriko)
set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(UpdateCacheVariable)
include(libigl)

# Work with non-standard homebrew installations
# (from ceres build system)
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

# ===  Try to find suitesparse === #
option(SUITESPARSE "Enable SuiteSparse." ON)
if (SUITESPARSE)
    message("-- Searching for SuiteSparse.")
    find_package(SuiteSparse)
    if (SUITESPARSE_FOUND)
        message("-- Found SuiteSparse ${SUITESPARSE_VERSION}, building with SuiteSparse.")
    else (SUITESPARSE_FOUND)
        message("-- Did not find all SuiteSparse dependencies, disabling SuiteSparse support.")
        update_cache_variable(SUITESPARSE OFF)
    endif (SUITESPARSE_FOUND)
else (SUITESPARSE)
    message("-- Building without SuiteSparse.")
endif (SUITESPARSE)

# HACK BY NICK
# The SuiteSparse logic above doesn't look for UMFpack, but we need it.
# This code attempts to find if by assuming that it will be in the same place
# as cholmod
if(SUITESPARSE AND SUITESPARSE_FOUND)
    string(REGEX REPLACE "cholmod" "umfpack" UMFPACK_LIBRARY ${CHOLMOD_LIBRARY})
    message("-- Guesstimated umfpack location as: ${UMFPACK_LIBRARY}")
    if(EXISTS ${UMFPACK_LIBRARY})
        list(APPEND SUITESPARSE_LIBRARIES ${UMFPACK_LIBRARY})
    else()
        message(WARNING "UMFPack guess failed, so we don't actually have SUITESPARSE support.")
        set(SUITESPARSE_FOUND FALSE)
    endif()
endif()

if(SUITESPARSE AND SUITESPARSE_FOUND)
    SET(GC_HAVE_SUITESPARSE TRUE)
    SET(GC_HAVE_SUITESPARSE TRUE)
    set(SUITESPARSE_INCLUDE_DIRS "${SUITESPARSE_INCLUDE_DIRS}")
    set(SUITESPARSE_LIBRARIES    "${SUITESPARSE_LIBRARIES}"   )
else()
    SET(GC_HAVE_SUITESPARSE FALSE)
    SET(GC_HAVE_SUITESPARSE FALSE PARENT_SCOPE)
endif()

# ===  done finding Suitesparse === #



add_library(${PROJECT_NAME} SHARED
        include/metriko/core/hmesh/hmesh.cpp
        include/metriko/core/vectorfield/face_rosy_field.cpp
        include/metriko/core/igm/parameterization_setup.cpp
        include/metriko/core/igm/parameterization_integ.cpp
        include/metriko/core/tmesh/motorcycle.cpp
        include/metriko/core/tmesh/tmesh.cpp
)
target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX _USE_MATH_DEFINES)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)
target_include_directories(${PROJECT_NAME} PUBLIC include)

if(GC_HAVE_SUITESPARSE)
    target_compile_definitions(${PROJECT_NAME} INTERFACE GC_HAVE_SUITESPARSE)
    target_include_directories(${PROJECT_NAME} INTERFACE ${SUITESPARSE_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} INTERFACE ${SUITESPARSE_LIBRARIES})
endif()
