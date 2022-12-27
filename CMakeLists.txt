#--------------------------------------------------------------
# CMake file for HydroChrono library and its demos
# Eventually will add SWIG and python lines to generate python
# library from C++ code
#--------------------------------------------------------------

cmake_minimum_required(VERSION 3.18.2)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# Default build type to use if none was specified
if(NOT DEFINED HYDRO_CHRONO_DEFAULT_BUILD_TYPE)
	set(HYDRO_CHRONO_DEFAULT_BUILD_TYPE "Release")
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to '${HYDRO_CHRONO_DEFAULT_BUILD_TYPE}' as none was specified.")

	set(CMAKE_BUILD_TYPE ${HYDRO_CHRONO_DEFAULT_BUILD_TYPE}
	    CACHE STRING "Choose the type of build." FORCE)
	mark_as_advanced(CMAKE_BUILD_TYPE)

	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#--------------------------------------------------------------
# Modify the project name if you want: 
#--------------------------------------------------------------

project(HydroChrono 
	VERSION 0.0.1
	DESCRIPTION "Hydro dynamic library based on chrono engine"
	LANGUAGES CXX
)


option (HYDRO_CHRONO_ENABLE_TESTS "Enable tests" ON)


# find required packages and libraries to make HydroChrono library
set (LIB_TYPE STATIC) # or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

find_package(Chrono COMPONENTS Irrlicht CONFIG)
find_package(HDF5 NAMES hdf5 COMPONENTS CXX ${SEARCH_TYPE})
#find_package(SWIG REQUIRED)
#find_package(PythonLibs)
#include(${SWIG_USE_FILE})

include_directories(${HDF5_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${CHRONO_INCLUDE_DIRS})
#include_directories(${HDF5_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${CHRONO_INCLUDE_DIRS})

set(PROJECT_LIBS ${HDF5_LIBS} ${HDF5_LIBRARIES})
set (LINK_LIBS ${LINK_LIBS} ${CHRONO_LIBRARIES} ${HDF5_CXX_${LIB_TYPE}_LIBRARY})
#SET(CMAKE_SWIG_FLAGS "")

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

#--------------------------------------------------------------
# === 3 ===
# Add the executable from your project and specify all C++ 
# files in your project. 
#--------------------------------------------------------------
set(HYDRO_CHRONO_SOURCES
  
	src/h5fileinfo.cpp
	src/chloadaddedmass.cpp
	src/hydro_forces.cpp

)

add_library(HydroChrono STATIC ${HYDRO_CHRONO_SOURCES} )


target_include_directories(
	HydroChrono

	PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>

	PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/
	${CMAKE_CURRENT_SOURCE_DIR}/include

	)




# add_executable(sphere_reg_waves "demos/sphere/sphere_reg_waves.cpp")
# add_executable(rm3_decay "demos/rm3/rm3_decay.cpp")
# add_executable(rm3_reg_waves "demos/rm3/rm3_reg_waves.cpp")
# add_executable(oswec_decay "demos/oswec/oswec_decay.cpp")
#SWIG_ADD_LIBRARY(py_wec_chrono LANGUAGE python SOURCES hydro_forces_no_viz.i)

target_compile_features(HydroChrono PUBLIC cxx_std_17)

#--------------------------------------------------------------
# Set properties for your executable target
# 
# Note that here we define a macro CHRONO_DATA_DIR which will
# contain the path to the Chrono data directory, either in its
# source tree (if using a build version of Chrono), or in its
# install tree (if using an installed version of Chrono).
#--------------------------------------------------------------

set_target_properties(HydroChrono PROPERTIES 
	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

# set_target_properties(sphere_decay PROPERTIES 
# 	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
# 	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
# 	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

# set_target_properties(sphere_reg_waves PROPERTIES 
# 	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
# 	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
# 	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

# set_target_properties(rm3_decay PROPERTIES 
# 	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
# 	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
# 	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

# set_target_properties(rm3_reg_waves PROPERTIES 
# 	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
# 	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
# 	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

# set_target_properties(oswec_decay PROPERTIES 
# 	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
# 	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
# 	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

#--------------------------------------------------------------
# Link to Chrono libraries and dependency libraries
#--------------------------------------------------------------

target_link_libraries(HydroChrono ${LINK_LIBS})

# target_link_libraries(sphere_reg_waves HydroChrono)
# target_link_libraries(rm3_decay HydroChrono)
# target_link_libraries(rm3_reg_waves HydroChrono)
# target_link_libraries(oswec_decay HydroChrono)
#SWIG_LINK_LIBRARIES(py_wec_chrono wec_chrono ${PYTHON_LIBRARIES})








add_executable(sphere_decay demos/sphere/sphere_decay.cpp)

target_link_libraries(sphere_decay HydroChrono)


# ====================
# TESTS
# ====================
if(HYDRO_CHRONO_ENABLE_TESTS)

	if(DEFINED ENV{HYDRO_CHRONO_DATA_DIR})
		# Override DATA_DIR by environement variable
		set(HYDRO_CHRONO_DATA_DIR $ENV{HYDRO_CHRONO_DATA_DIR})
	else()
		set(HYDRO_CHRONO_DATA_DIR "${PROJECT_SOURCE_DIR}/demos")
	endif()

	include(CTest)
	enable_testing()

	add_subdirectory(tests)

endif(HYDRO_CHRONO_ENABLE_TESTS)


#--------------------------------------------------------------
# === 4 (OPTIONAL) ===
# 
# Optionally, add a custom command for copying all Chrono and
# dependency DLLs to the appropriate binary output folder.
# This function has effect only on Windows.
# 
# Note that you must first set EXECUTABLE_OUTPUT_PATH
# (this can simply be ${CMAKE_BINARY_DIR}, like in this example)
#--------------------------------------------------------------

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#add_DLL_copy_command("${CHRONO_DLLS}")
