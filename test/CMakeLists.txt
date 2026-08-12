cmake_minimum_required(VERSION 3.28)
project(metriko_test)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# metriko を自分で取り込む独立プロジェクト（example/ と同じ形）
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. "metriko")
enable_testing()

add_executable(test_tmesh test_tmesh.cpp)
target_link_libraries(test_tmesh PRIVATE metriko)
target_include_directories(test_tmesh PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# 全 tquad の region(allowed_range) を JSON にダンプするツール
add_executable(dump_region dump_region.cpp)
target_link_libraries(dump_region PRIVATE metriko)
target_include_directories(dump_region PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# allowed_range regression checker: compares against test/fixtures/tquad_region/<mesh>.json
add_executable(test_allowed_range test_allowed_range.cpp)
target_link_libraries(test_allowed_range PRIVATE metriko)
target_include_directories(test_allowed_range PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# test settings come from test_config.json (read at configure time).
#   { "directory": "/abs/dir",
#     "patterns": [ { "gridscale":0.05, "vectorfield":"smoothest"|"curvature_aligned",
#                     "expected":"expected_allowed.json", "meshes":["torus","..."] }, ... ] }
#   mesh path = directory + "/" + <name> + ".obj"; tests registered per (vectorfield, mesh).
#   region fixture (optional) at fixtures/tquad_region/<vectorfield>/<name>.json.
#   * editing the JSON requires a CMake reconfigure (auto on cmake --build / CLion reload).
set(METRIKO_TEST_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/test_config.json"
    CACHE FILEPATH "test settings json (directory, patterns[gridscale,vectorfield,expected,meshes])")
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
                add_test(NAME tmesh.${_vf}.${_name} COMMAND test_tmesh ${_gs} ${_vf} ${_mesh})
                set(_fx "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/tquad_region/${_vf}/${_name}.json")
                if(EXISTS ${_fx})
                    add_test(NAME region.${_vf}.${_name} COMMAND test_allowed_range ${_gs} ${_vf} ${_mesh} ${_fx})
                endif()
            endforeach()
        endforeach()
    else()
        message(WARNING "test config '${METRIKO_TEST_CONFIG}' has empty 'patterns'; no tests registered.")
    endif()
endif()
