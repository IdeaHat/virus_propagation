cmake_minimum_required(VERSION 2.8)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_subdirectory(graph_manipulation)

PROJECT(calc_sis)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)
find_package(Eigen3 REQUIRED)
include_directories(${graph_manipulation_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})
ADD_EXECUTABLE(${PROJECT_NAME} main.cpp)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} graph_manipulation)
