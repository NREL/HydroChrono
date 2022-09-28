
cmake_minimum_required(VERSION 3.18.2)
project(demos CXX)

add_executable(sphere_decay
  ${CMAKE_SOURCE_DIR}/demos/sphere/sphere_decay.cpp
  ${CMAKE_SOURCE_DIR}/src/hydro_forces.cpp
)

target_include_directories(sphere_decay PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/include
 )