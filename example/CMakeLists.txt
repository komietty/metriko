cmake_minimum_required(VERSION 3.28)
project(example)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. "metriko")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/polyscope/ "polyscope")
add_executable(${PROJECT_NAME} example_1.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE metriko polyscope)
