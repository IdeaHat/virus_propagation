cmake_minimum_required(VERSION 2.8)
PROJECT(graph_manipulation)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)
find_package(Eigen3 REQUIRED)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}" "${EIGEN3_INCLUDE_DIR}")
ADD_LIBRARY(${PROJECT_NAME} 
  "graph_manip/graph_data_structures.hpp"
  graph_data_structures.cpp)

