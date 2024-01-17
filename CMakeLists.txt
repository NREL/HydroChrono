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
if(NOT DEFINED HYDROCHRONO_DEFAULT_BUILD_TYPE)
	set(HYDROCHRONO_DEFAULT_BUILD_TYPE "Release")
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to '${HYDROCHRONO_DEFAULT_BUILD_TYPE}' as none was specified.")

	set(CMAKE_BUILD_TYPE ${HYDROCHRONO_DEFAULT_BUILD_TYPE}
	    CACHE STRING "Choose the type of build." FORCE)
	mark_as_advanced(CMAKE_BUILD_TYPE)

	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#--------------------------------------------------------------
# Modify the project name if you want: 
#--------------------------------------------------------------

project(HydroChrono 
	VERSION 0.1
	DESCRIPTION "A hydrodynamic library for use with Project Chrono."
	LANGUAGES CXX
)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# ===================
# OPTIONS
# ===================
option (HYDROCHRONO_ENABLE_TESTS "Enable tests" ON)
option (HYDROCHRONO_ENABLE_IRRLICHT "Enable irrlicht visualization library" ON)
option (HYDROCHRONO_ENABLE_DEMOS "Enable demo executables" ON)
option (HYDROCHRONO_ENABLE_USER_DOC "User's documentation" OFF)
option (HYDROCHRONO_ENABLE_PROG_DOC "Programmer's documentation" OFF)


# find required packages and libraries to make HydroChrono library
set (LIB_TYPE STATIC) # or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	# If irrlicht is required chrono library must be compiled with irrlicht
	find_package(Chrono COMPONENTS Irrlicht CONFIG REQUIRED)
	add_compile_definitions(HYDROCHRONO_HAVE_IRRLICHT=1)
else(HYDROCHRONO_ENABLE_IRRLICHT)
	find_package(Chrono CONFIG REQUIRED)
endif(HYDROCHRONO_ENABLE_IRRLICHT)

find_package(HDF5 NAMES hdf5 COMPONENTS CXX ${SEARCH_TYPE})


#-----------------------------------------------------------------------------
# Fix for VS 2017 15.8 and newer to handle alignment specification with Eigen
#-----------------------------------------------------------------------------

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  if(MSVC AND ${MSVC_VERSION} GREATER_EQUAL 1915)
    add_compile_definitions( "ENABLE_EXTENDED_ALIGNED_STORAGE" )
  endif()
endif()

#--------------------------------------------------------------
# Tweaks to disable some warnings with MSVC
#--------------------------------------------------------------
if(MSVC)
	# avoids deprecation warnings
    add_definitions(
		-D_CRT_SECURE_NO_DEPRECATE  
		-D_SCL_SECURE_NO_DEPRECATE
		-DNOMINMAX
	)
    add_compile_options(/wd4275)                   # disable warnings triggered by Irrlicht
    add_compile_options(/wd4251)                   # class needs to have dll-interface
endif()

# =======================
# HydroChrono Library
# =======================
set(HYDROCHRONO_SOURCES
  
	src/h5fileinfo.cpp
	src/chloadaddedmass.cpp
	src/hydro_forces.cpp
	src/helper.cpp
	src/wave_types.cpp

)

add_library(HydroChrono STATIC ${HYDROCHRONO_SOURCES} )

target_compile_features(HydroChrono PUBLIC cxx_std_17)

target_include_directories(
	HydroChrono

	PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
	${CHRONO_INCLUDE_DIRS} 	

	PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/
	${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_definitions(HydroChrono

	PUBLIC

	PRIVATE 
		CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"
)

target_compile_options(HydroChrono BEFORE 
	PUBLIC
		${CHRONO_CXX_FLAGS}
)

target_link_options(HydroChrono BEFORE 
	PUBLIC
		${CHRONO_LINKER_FLAGS}
)

set_target_properties(HydroChrono
	PROPERTIES
	POSITION_INDEPENDENT_CODE ON
)

# temporary hack
# seems there is no namespace target hdf5:: on linux
# for HDF5 version 1.10 to 1.14 ?
if(NOT TARGET hdf5::hdf5_cpp-static)
	add_library(hdf5::hdf5_cpp-static ALIAS hdf5_cpp-static)
endif()

target_link_libraries(HydroChrono 
	PUBLIC
		${CHRONO_LIBRARIES}		
	PRIVATE
		hdf5::hdf5_cpp-static

)

# ====================
# Irrlicht GUI helper
# ====================

add_library(HydroChronoGUI)

target_sources(
	HydroChronoGUI

	PUBLIC
	src/gui/guihelper.cpp
)

target_compile_features(HydroChronoGUI PUBLIC cxx_std_17)


target_include_directories(
	HydroChronoGUI

	PUBLIC
		$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:include>
		${CHRONO_INCLUDE_DIRS}
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
endif(HYDROCHRONO_ENABLE_IRRLICHT)

target_compile_options(HydroChronoGUI BEFORE
	PUBLIC
		${CHRONO_CXX_FLAGS}
)

target_link_options(HydroChronoGUI BEFORE
	PUBLIC
		${CHRONO_LINKER_FLAGS}
)

target_link_libraries(HydroChronoGUI
	PRIVATE
		${CHRONO_LIBRARIES}
)



# ====================
# DEMOS
# ====================
if(HYDROCHRONO_ENABLE_DEMOS)
	add_subdirectory(demos)
endif(HYDROCHRONO_ENABLE_DEMOS)


# ====================
# TESTS
# ====================
if(HYDROCHRONO_ENABLE_TESTS)

	if(DEFINED ENV{HYDROCHRONO_DATA_DIR})
		# Override DATA_DIR by environement variable
		set(HYDROCHRONO_DATA_DIR $ENV{HYDROCHRONO_DATA_DIR})
	else()
		set(HYDROCHRONO_DATA_DIR "${PROJECT_SOURCE_DIR}/demos")
	endif()

	include(CTest)
	enable_testing()

	add_subdirectory(tests)

endif(HYDROCHRONO_ENABLE_TESTS)



# Not a good idea to copy DLLs
# dosen't work on Linux
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#add_DLL_copy_command("${CHRONO_DLLS}")

include(GNUInstallDirs)

# Export to include as external library in other projects

install(TARGETS HydroChrono
    EXPORT HydroChronoTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    EXPORT HydroChronoTargets
    FILE HydroChronoTargets.cmake
    DESTINATION "${CMAKE_INSTALL_DATADIR}/HydroChrono/cmake"
    NAMESPACE HydroChrono::
)

install(
    DIRECTORY include
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/HydroChronoConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/HydroChronoConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_DATADIR}/HydroChrono/cmake"
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/HydroChronoConfig.cmake"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/HydroChrono/cmake"
)

export(
    EXPORT HydroChronoTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/HydroChronoTargets.cmake"
    NAMESPACE HydroChrono::
)
