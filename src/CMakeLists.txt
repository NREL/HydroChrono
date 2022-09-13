# CMakeLists.txt in src folder

option(BUILD_DEMOS "Build demo programs" ON)

find_package(Chrono COMPONENTS Irrlicht CONFIG)

find_package(HDF5 NAMES hdf5 COMPONENTS CXX ${SEARCH_TYPE})
    if(NOT HDF5_FOUND)
        # Look for a package configuration file
        # LIB_TYPE can be one of STATIC or SHARED.
        set(LIB_TYPE SHARED) # STATIC or SHARED
        string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)
        find_package(HDF5 NAMES hdf5 COMPONENTS CXX ${SEARCH_TYPE})

        if (HDF5_FOUND)
            set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
            set(HDF5_C_LIBRARIES ${HDF5_C_${LIB_TYPE}_LIBRARY})
            set(HDF5_CXX_LIBRARIES ${HDF5_CXX_${LIB_TYPE}_LIBRARY)
        endif()
    endif()


#find_package(SWIG REQUIRED)
#find_package(PythonLibs)
#include(${SWIG_USE_FILE})
# find required packages and libraries to make HydroChrono library
# set (LIB_TYPE STATIC) # or SHARED
# string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

set(PROJECT_LIBS ${HDF5_LIBS} ${HDF5_LIBRARIES})
set (LINK_LIBS ${LINK_LIBS} ${CHRONO_LIBRARIES} ${HDF5_CXX_${LIB_TYPE}_LIBRARY})
#SET(CMAKE_SWIG_FLAGS "")

#-----------------------------------------------------------------------------
# Set the base include directories
#-----------------------------------------------------------------------------

include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CHRONO_INCLUDE_DIRS})
#include_directories(${PYTHON_INCLUDE_PATH})

#--------------------------------------------------------------
# === 3 ===
# Add the executable from your project and specify all C++ 
# files in your project. 
#--------------------------------------------------------------

add_library(HydroChrono STATIC "hydro_forces.cpp" "hydro_forces.h")
add_executable(sphere_decay "../demos/sphere/sphere_decay.cpp")
add_executable(sphere_reg_waves "../demos/sphere/sphere_reg_waves.cpp")
add_executable(rm3_decay "../demos/rm3/rm3_decay.cpp")
add_executable(rm3_reg_waves "../demos/rm3/rm3_reg_waves.cpp")
#SWIG_ADD_LIBRARY(py_wec_chrono LANGUAGE python SOURCES hydro_forces_no_viz.i)


target_compile_features(HydroChrono PUBLIC cxx_std_17)

#--------------------------------------------------------------
# Link to Chrono libraries and dependency libraries
#--------------------------------------------------------------

target_link_libraries(HydroChrono ${LINK_LIBS})
target_link_libraries(sphere_decay HydroChrono)
target_link_libraries(sphere_reg_waves HydroChrono)
target_link_libraries(rm3_decay HydroChrono)
target_link_libraries(rm3_reg_waves HydroChrono)
#SWIG_LINK_LIBRARIES(py_wec_chrono wec_chrono ${PYTHON_LIBRARIES})


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

set_target_properties(sphere_decay PROPERTIES 
	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

set_target_properties(sphere_reg_waves PROPERTIES 
	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

set_target_properties(rm3_decay PROPERTIES 
	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")

set_target_properties(rm3_reg_waves PROPERTIES 
	    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}"
	    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
	    LINK_FLAGS "${CHRONO_LINKER_FLAGS}")



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

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
add_DLL_copy_command("${CHRONO_DLLS}")
# if(BUILD_DEMOS)
#   add_subdirectory(../demos)
# endif()