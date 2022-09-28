
cmake_minimum_required(VERSION 3.18.2)
project(demos CXX)

#######################################################
# demo/sphere
#######################################################
# add_executable for sphere demo, list it's source files and ${CMAKE_SOURCE_DIR}/src/hydro_forces.cpp
add_executable(sphere_decay
  ${CMAKE_SOURCE_DIR}/demos/sphere/sphere_decay.cpp
  ${CMAKE_SOURCE_DIR}/src/hydro_forces.cpp
)

# include directories for shere demo 
target_include_directories(sphere_decay PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/include
 )

# link hydro_lib and dependent libraries to sphere_decay demo 
target_link_libraries(sphere_decay PUBLIC hydro_lib)

# set some options for chrono demos here
set_target_properties(sphere_decay PROPERTIES 
    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

#######################################################
# demo/newdemo
# example lines to add your own demo to the cmake 
#######################################################
# # add_executable for newdemo demo, list it's source files and ${CMAKE_SOURCE_DIR}/src/hydro_forces.cpp
# add_executable(newdemo
#   ${CMAKE_SOURCE_DIR}/demos/newdemo/newdemo.cpp
#   ${CMAKE_SOURCE_DIR}/src/hydro_forces.cpp
# )

# # include directories for newdemo demo 
# target_include_directories(newdemo PUBLIC
#   ${CMAKE_SOURCE_DIR}/src
#   ${CMAKE_SOURCE_DIR}/include
#  )

# # link hydro_lib and dependent libraries to newdemo demo 
# target_link_libraries(newdemo PUBLIC hydro_lib)

# # set some options for chrono demos here
# set_target_properties(newdemo PROPERTIES 
#     COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
#     COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
#     LINK_FLAGS "${CHRONO_LINKER_FLAGS}")