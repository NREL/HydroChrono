# CMakeLists.txt in src folder

cmake_minimum_required(VERSION 3.18.2)
project(target-system LANGUAGES CXX)

find_package(Chrono COMPONENTS Irrlicht CONFIG)

set(LIB_TYPE STATIC) # STATIC or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)
find_package(HDF5 NAMES hdf5 COMPONENTS CXX ${SEARCH_TYPE})
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}")
set (LINK_LIBS ${LINK_LIBS} ${HDF5_CXX_${LIB_TYPE}_LIBRARY})
    # if(NOT HDF5_FOUND)
    #     # Look for a package configuration file
    #     # LIB_TYPE can be one of STATIC or SHARED.
    #     set(LIB_TYPE STATIC) # STATIC or SHARED
    #     string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)
    #     find_package(HDF5 NAMES hdf5 COMPONENTS CXX ${SEARCH_TYPE})
    #     if (HDF5_FOUND)
    #         set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR})
    #         set(HDF5_C_LIBRARIES ${HDF5_C_${LIB_TYPE}_LIBRARY})
    #         set(HDF5_CXX_LIBRARIES ${HDF5_CXX_${LIB_TYPE}_LIBRARY})
    #     endif()
    # endif()
#list cpp files for library here
add_library(system STATIC # static here? or object?
	${CMAKE_SOURCE_DIR}/src/hydro_forces.cpp
)
#list directories with .h files here
target_include_directories(system PUBLIC
	${CMAKE_SOURCE_DIR}/include
)
#include external dependencies
target_include_directories(system PUBLIC 
	${HDF5_INCLUDE_DIRS}
	${CHRONO_INCLUDE_DIRS}
)

# this fixes unresolved chrono link stuff, but not hdf5 stuff
message("HDF5_CXX_${LIB_TYPE}_LIBRARY = ${HDF5_CXX_${LIB_TYPE}_LIBRARY}")
set (LINK_LIBS ${LINK_LIBS} ${CHRONO_LIBRARIES} ${HDF5_CXX_${LIB_TYPE}_LIBRARY})
target_link_libraries(system ${LINK_LIBS})



#find_package(SWIG REQUIRED)
#find_package(PythonLibs)
#include(${SWIG_USE_FILE})
# find required packages and libraries to make HydroChrono library
# set (LIB_TYPE STATIC) # or SHARED
# string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

#SET(CMAKE_SWIG_FLAGS "")


#include_directories(${PYTHON_INCLUDE_PATH})

#SWIG_ADD_LIBRARY(py_wec_chrono LANGUAGE python SOURCES hydro_forces_no_viz.i)
#SWIG_LINK_LIBRARIES(py_wec_chrono wec_chrono ${PYTHON_LIBRARIES})
