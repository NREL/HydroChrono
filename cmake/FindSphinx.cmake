# From https://www.vortech.nl/en/integrating-sphinx-in-cmake/

include(FindPackageHandleStandardArgs)

# We are likely to find Sphinx near the Python interpreter
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    get_filename_component(_PYTHON_DIR "${Python3_EXECUTABLE}" DIRECTORY)
    set(
        _PYTHON_PATHS
        "${_PYTHON_DIR}"
        "${_PYTHON_DIR}/bin"
        "${_PYTHON_DIR}/Scripts")
endif()

find_program(
    SPHINX_EXECUTABLE
    NAMES sphinx-build sphinx-build.exe
    HINTS ${_PYTHON_PATHS})
mark_as_advanced(SPHINX_EXECUTABLE)

find_package_handle_standard_args(Sphinx DEFAULT_MSG SPHINX_EXECUTABLE)

# If finding Sphinx fails, there is no use in defining
# add_sphinx_document, so return early
if(NOT Sphinx_FOUND)
    return()
endif()

# add_sphinx_document(
#   <name>
#   CONF_FILE <conf-py-filename>
#   [C_API <c-api-header-file>]
#   [SKIP_HTML] [SKIP_PDF]
#   <rst-src-file>...)
#
# Function for creating Sphinx documentation targets.
function(add_sphinx_document TARGET_NAME)

    cmake_parse_arguments(
        ${TARGET_NAME}
        "SKIP_HTML;SKIP_PDF"
        "CONF_FILE"
        ""
        ${ARGN})

        get_filename_component(SRCDIR "${${TARGET_NAME}_CONF_FILE}" DIRECTORY)
        set(INTDIR "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/source")
        set(OUTDIR "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/build")

        string(TIMESTAMP SPHINX_TARGET_YEAR "%Y" UTC)

        add_custom_command(
            OUTPUT "${INTDIR}/conf.py"
            COMMAND "${CMAKE_COMMAND}" -E make_directory "${INTDIR}"
            COMMAND
                "${CMAKE_COMMAND}"
                "-DCONFIGURE_FILE_IN=${${TARGET_NAME}_CONF_FILE}"
                "-DCONFIGURE_FILE_OUT=${INTDIR}/conf.py"
                "-DSPHINX_TARGET_NAME=${TARGET_NAME}"
                "-DSPHINX_TARGET_VERSION=${PROJECT_VERSION}"
                "-DSPHINX_TARGET_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}"
                "-DSPHINX_TARGET_VERSION_MINOR=${PROJECT_VERSION_MINOR}"
                "-DSPHINX_TARGET_YEAR=${SPHINX_TARGET_YEAR}"
                -P "${_SPHINX_SCRIPT_DIR}/BuildTimeConfigureFile.cmake"
            DEPENDS "${${TARGET_NAME}_CONF_FILE}")
        
        set(SPHINX_DEPENDS "${INTDIR}/conf.py") 
        
        set(_SPHINX_SCRIPT_DIR ${CMAKE_CURRENT_LIST_DIR})

endfunction()