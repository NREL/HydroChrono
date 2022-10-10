#--------------------------------------------------------------
# CMake file for HydroChrono library and its demos
# Eventually will add SWIG and python lines to generate python
# library from C++ code
#--------------------------------------------------------------

cmake_minimum_required(VERSION 3.18.2)

# set_property(GLOBAL PROPERTY USE_FOLDERS Off)

project(HydroChrono LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

#-----------------------------------------------------------------------------
# Fix for VS 2017 15.8 and newer to handle alignment specification with Eigen
#-----------------------------------------------------------------------------
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  if(MSVC AND ${MSVC_VERSION} GREATER_EQUAL 1915)
    add_definitions( "-D_ENABLE_EXTENDED_ALIGNED_STORAGE" )
  endif()
endif()

#--------------------------------------------------------------
# Tweaks to disable some warnings with MSVC
#--------------------------------------------------------------
if(MSVC)
    add_definitions("-D_CRT_SECURE_NO_DEPRECATE")  # avoids deprecation warnings
    add_definitions("-D_SCL_SECURE_NO_DEPRECATE")  # avoids deprecation warnings
    add_definitions( "-DNOMINMAX" )                # do not use MSVC's min/max macros
    add_compile_options(/wd4275)                   # disable warnings triggered by Irrlicht
    add_compile_options(/wd4251)                   # class needs to have dll-interface
else()
    set(EXTRA_COMPILE_FLAGS "")
endif()

#-----------------------------------------------------------------------------
# Specify locations of subdirectories
#-----------------------------------------------------------------------------
add_subdirectory(${CMAKE_SOURCE_DIR}/src)


option(BUILD_DEMOS "Build demo programs" ON)
# if option for demos on:
if(BUILD_DEMOS AND IS_DIRECTORY ${CMAKE_SOURCE_DIR}/demos)
	add_subdirectory(${CMAKE_SOURCE_DIR}/demos)

# # copy input files
# add_custom_command(
#         TARGET sphere_decay POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy
#                 ${CMAKE_SOURCE_DIR}/demos/sphere/hydroData/sphere.h5
#                 ${CMAKE_CURRENT_BINARY_DIR}/demos/sphere.h5)

endif()