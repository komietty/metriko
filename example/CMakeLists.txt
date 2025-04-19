cmake_minimum_required(VERSION 3.28)
project(metriko_example)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. "metriko")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../libs/polyscope/ "polyscope")
add_executable(${PROJECT_NAME} example_qgp.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE metriko polyscope)
