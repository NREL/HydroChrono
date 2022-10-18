project(target-hydropy LANGUAGES CXX)

find_package(SWIG REQUIRED) # this adds FindSWIG cmake stuff
include(${SWIG_USE_FILE}) # this loads cmake file swig_use_file  which adds swig_add_library
message("SWIG_USE_FILE = ${SWIG_USE_FILE}")
if(NOT SWIG_FOUND)
	message("SWIG NOT FOUND")
endif()

# find_package(Python)
# set(Python_ADDITIONAL_VERSIONS 3.10.4)
find_package(Python COMPONENTS Development)

message("Python_FOUND:${Python_FOUND}")
message("Python_VERSION:${Python_VERSION}")
message("Python_Development_FOUND:${Python_Development_FOUND}")
message("Python_LIBRARIES:${Python_LIBRARIES}")

include_directories(${Python_INCLUDE_DIRS}) # https://cmake.org/cmake/help/v3.12/module/FindPython.html

message("Python_INCLUDE_DIRS = ${Python_INCLUDE_DIRS}")
message("Python_LIBRARIES = ${Python_LIBRARIES}")
if(NOT ${Python_FOUND})
	message("Python NOT FOUND")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include
					${CMAKE_SOURCE_DIR}/src/hydro_cxx	
					${CMAKE_SOURCE_DIR}/src/hydro_swig)

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(swig_test.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(swig_test.i PROPERTIES SWIG_FLAGS "-includeall")

swig_add_library(hydro_py TYPE SHARED LANGUAGE python SOURCES swig_test.i)
target_include_directories(hydro_py PRIVATE ${Python_INCLUDE_DIRS})
swig_link_libraries(hydro_py ${Python_LIBRARIES} ${LINK_LIBS} hydro_lib)

set_target_properties(${SWIG_MODULE_hydro_py_REAL_NAME} PROPERTIES
                      PROJECT_LABEL "HydroChrono_python_hydro_py"
                      OUTPUT_NAME   "${SWIG_MODULE_hydro_py_REAL_NAME}"
                      LIBRARY_OUTPUT_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
                      )