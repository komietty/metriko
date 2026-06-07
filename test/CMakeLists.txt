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

# テスト対象メッシュは構成時にカンマ区切り・絶対パスで指定:
#   cmake -DMETRIKO_TEST_MESHES=/abs/a.obj,/abs/b.obj
# CLion の CMake options 欄でも同様（カンマ区切り。';'/'\;' のエスケープ不要）。
set(METRIKO_TEST_MESHES    ""     CACHE STRING "mesh paths to test (comma-separated, absolute)")
set(METRIKO_TEST_GRIDSCALE "0.05" CACHE STRING "grid scale for tmesh test")

if(NOT METRIKO_TEST_MESHES)
    message(WARNING "METRIKO_TEST_MESHES is empty; no mesh tests registered. "
                    "Set e.g. -DMETRIKO_TEST_MESHES=/abs/a.obj,/abs/b.obj")
endif()

string(REPLACE "," ";" _meshes "${METRIKO_TEST_MESHES}")
foreach(mesh ${_meshes})
    get_filename_component(name ${mesh} NAME_WE)
    add_test(NAME tmesh.${name} COMMAND test_tmesh ${METRIKO_TEST_GRIDSCALE} ${mesh})
endforeach()
