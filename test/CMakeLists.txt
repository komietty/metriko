cmake_minimum_required(VERSION 3.28)
project(metriko_test)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "-O2 -g")

# metriko を自分で取り込む独立プロジェクト（example/ と同じ形）
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. "metriko")
enable_testing()

# end-to-end smoke test through quad extraction (lib.h)
add_executable(test_remesh test_remesh.cpp)
target_link_libraries(test_remesh PRIVATE metriko)
target_include_directories(test_remesh PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# test settings come from test_config.json (read at configure time).
#   { "directory": "/abs/dir",
#     "patterns": [ { "gridscale":0.05, "vectorfield":"smoothest"|"curvature_aligned",
#                     "meshes":["torus","..."] }, ... ] }
#   mesh path = directory + "/" + <name> + ".obj"; tests registered per (vectorfield, mesh).
#   * editing the JSON requires a CMake reconfigure (auto on cmake --build / CLion reload).
set(METRIKO_TEST_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/test_config.json"
    CACHE FILEPATH "test settings json (directory, patterns[gridscale,vectorfield,meshes])")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${METRIKO_TEST_CONFIG})

if(NOT EXISTS ${METRIKO_TEST_CONFIG})
    message(WARNING "test config not found: ${METRIKO_TEST_CONFIG}; no tests registered.")
else()
    file(READ ${METRIKO_TEST_CONFIG} _cfg)
    string(JSON _dir GET ${_cfg} directory)
    string(JSON _np  LENGTH ${_cfg} patterns)
    if(_np GREATER 0)
        math(EXPR _lastp "${_np} - 1")
        foreach(_pi RANGE ${_lastp})
            string(JSON _p  GET ${_cfg} patterns ${_pi})
            string(JSON _gs GET ${_p} gridscale)
            string(JSON _vf GET ${_p} vectorfield)
            string(JSON _nm LENGTH ${_p} meshes)
            if(NOT _nm GREATER 0)
                continue()
            endif()
            math(EXPR _lastm "${_nm} - 1")
            foreach(_mi RANGE ${_lastm})
                string(JSON _name GET ${_p} meshes ${_mi})
                set(_mesh "${_dir}/${_name}.obj")
                # test names carry the gridscale: the same mesh may appear in several patterns
                set(_id "${_vf}.${_name}.${_gs}")
                add_test(NAME remesh.${_id} COMMAND test_remesh ${_gs} ${_vf} ${_mesh})
                # schedule heavy meshes first: file size is a good enough proxy for run time
                if(EXISTS ${_mesh})
                    file(SIZE ${_mesh} _sz)
                    set_tests_properties(remesh.${_id} PROPERTIES COST ${_sz})
                endif()
            endforeach()
        endforeach()
    else()
        message(WARNING "test config '${METRIKO_TEST_CONFIG}' has empty 'patterns'; no tests registered.")
    endif()
endif()
