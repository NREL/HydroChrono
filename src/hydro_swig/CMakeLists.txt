project(target-hydropy LANGUAGES CXX)

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH}) # or ${PYTHON_INCLUDE_DIR}?

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/include
					${CMAKE_SOURCE_DIR}/src/hydro_cxx	
					${CMAKE_SOURCE_DIR}/src/hydro_swig)

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(swig_test.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(swig_test.i PROPERTIES SWIG_FLAGS "-includeall")

swig_add_library(hydro_py TYPE SHARED LANGUAGE python SOURCES swig_test.i)
target_include_directories(hydro_py PRIVATE ${PYTHON_INCLUDE_DIR})
swig_link_libraries(hydro_py ${PYTHON_LIBRARIES})