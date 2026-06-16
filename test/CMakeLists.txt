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

# テスト設定は test_config.json から CMake が構成時に読み、メッシュ毎に登録する。
#   test_config.json: { "gridscale": 0.05, "models_dir": "/abs/dir", "meshes": ["torus", ...] }
#   ・メッシュパスは models_dir + "/" + <name> + ".obj"。
#   ※ JSON を編集したら CMake 再構成が必要（cmake --build で自動 / CLion は test プロジェクトを Reload）。
set(METRIKO_TEST_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/test_config.json"
    CACHE FILEPATH "test settings json (gridscale, models_dir, meshes)")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${METRIKO_TEST_CONFIG})

if(NOT EXISTS ${METRIKO_TEST_CONFIG})
    message(WARNING "test config not found: ${METRIKO_TEST_CONFIG}; no tests registered.")
else()
    file(READ ${METRIKO_TEST_CONFIG} _cfg)
    string(JSON _gridscale  GET ${_cfg} gridscale)
    string(JSON _models_dir GET ${_cfg} models_dir)
    string(JSON _n_meshes LENGTH ${_cfg} meshes)
    if(_n_meshes GREATER 0)
        math(EXPR _last "${_n_meshes} - 1")
        foreach(_i RANGE ${_last})
            string(JSON _name GET ${_cfg} meshes ${_i})
            add_test(NAME tmesh.${_name} COMMAND test_tmesh ${_gridscale} ${_models_dir}/${_name}.obj)
        endforeach()
    else()
        message(WARNING "test config '${METRIKO_TEST_CONFIG}' has empty 'meshes'; no tests registered.")
    endif()

    # region regression: one test per saved fixture (mesh name = fixture basename).
    # gridscale must match the value the fixture was generated with (= test_config gridscale).
    file(GLOB _fixtures CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/tquad_region/*.json")
    foreach(_fx ${_fixtures})
        get_filename_component(_name ${_fx} NAME_WE)
        add_test(NAME region.${_name} COMMAND test_allowed_range ${_gridscale} ${_models_dir}/${_name}.obj ${_fx})
    endforeach()
endif()
