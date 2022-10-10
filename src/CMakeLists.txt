

add_subdirectory(${CMAKE_SOURCE_DIR}/src/hydro_cxx)

option(BUILD_PYTHON "Build Python API" ON)
# if option for python on:
if(BUILD_PYTHON AND IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src/hydro_swig)
	add_subdirectory(${CMAKE_SOURCE_DIR}/src/hydro_swig)
endif()