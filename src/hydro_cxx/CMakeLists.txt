# CMakeLists.txt in src folder

cmake_minimum_required(VERSION 3.18.2)
project(target-hydro_lib LANGUAGES CXX)

# find project chrono
find_package(Chrono COMPONENTS Irrlicht CONFIG)

# find hdf5 library
set(LIB_TYPE STATIC) # STATIC or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)
find_package(HDF5 NAMES hdf5 COMPONENTS CXX ${SEARCH_TYPE})
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}")

#list cpp files for library here
add_library(hydro_lib STATIC
	${CMAKE_SOURCE_DIR}/src/hydro_cxx/hydro_forces.cpp
	${CMAKE_SOURCE_DIR}/src/hydro_cxx/swig_test.cpp
)

#list directories with .h files here
target_include_directories(hydro_lib PUBLIC
	${CMAKE_SOURCE_DIR}/include
)

#include external dependencies like hdf5 and chrono
target_include_directories(hydro_lib PUBLIC 
	${HDF5_INCLUDE_DIRS}
	${CHRONO_INCLUDE_DIRS}
)

# link chrono and hdf5 to hydro_lib
set (LINK_LIBS ${LINK_LIBS} ${CHRONO_LIBRARIES} ${HDF5_CXX_${LIB_TYPE}_LIBRARY})
target_link_libraries(hydro_lib ${LINK_LIBS})

set_target_properties(hydro_lib PROPERTIES 
    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")


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
