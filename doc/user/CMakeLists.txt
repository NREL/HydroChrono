find_package(Python3 COMPONENTS Interpreter)

if(NOT Python3_FOUND)
    message(WARNING "Can't build user docs. Python3 not found!")
    message(WARNING "Disable user doc option or provide a Python3 environment") 
    message(WARNING "or provide hints in the CMakeCache.txt file: Python_ROOT_DIR=")       
    message(ERROR "and for a virtual env.  Python_FIND_VIRTUALENV=FIRST see cmake documentation")
else(NOT Python3_FOUND)
    message(STATUS "Found python: ${Python3_EXECUTABLE}")
endif(NOT Python3_FOUND)


find_package(Sphinx)
if(NOT Sphinx_FOUND)
    message(WARNING "Can't build user's docs. Sphinx not found!")
    message(ERROR "Disable user doc option or provide a sphinx environment")    
else(NOT Sphinx_FOUND)
    message(STATUS "Found sphinx: ${SPHINX_EXECUTABLE}")
endif(NOT Sphinx_FOUND)



# Check for additional modules necessary for sphinx documentation: breathe, ...
# Modules are configured in 'conf.py.in' sphinx configuration file
# Define variables for python code (booolean with integer)
set(PYHAVE_BREATHE 1)
set(PYHAVE_SPINXCONTRIB_BIBTEX 1)
if(HYDROCHRONO_ENABLE_PROG_DOC)
    set(PYHAVE_PROG_DOC 1)
else()
    set(PYHAVE_PROG_DOC 0)
endif()

if (Python3_FOUND AND Sphinx_FOUND)
    # mandatory bibtex contrib
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -E -c "import sphinxcontrib.bibtex"
        RESULT_VARIABLE HAVENOT_SPINXCONTRIB_BIBTEX
        OUTPUT_QUIET    
    )


    if (NOT ${HAVENOT_SPINXCONTRIB_BIBTEX} EQUAL 0)       
        set(PYHAVE_SPINXCONTRIB_BIBTEX  0)
        message(
                FATAL_ERROR
                "The \"sphinxcontrib.bibtex\" Python3 package is not installed. Please install it using the following command: \"pip3 install sphinxcontrib.bibtex\"."
        ) 
    else()
        message(STATUS "Found python 'sphinxcontrib.bibtex'")           
    endif()

    # optional breathe for doxygen
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -E -c "import breathe"
        RESULT_VARIABLE HAVENOT_BREATHE
        OUTPUT_QUIET    
    )
    if (NOT ${HAVENOT_BREATHE} EQUAL 0)
        set(PYHAVE_BREATHE  0)
    else()
        message(STATUS "Found python 'breathe'")
    endif()
endif()



set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/html)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

configure_file(conf.py.in ${SPHINX_BUILD}/conf.py @ONLY)


add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}

    COMMAND

    ${SPHINX_EXECUTABLE} -a -b html -c ${SPHINX_BUILD} # Force write all files
    ${SPHINX_SOURCE} ${SPHINX_BUILD}

    WORKING_DIRECTORY

    ${CMAKE_CURRENT_BINARY_DIR}

    DEPENDS

    ${SPHINX_SOURCE} # Every source files

    MAIN_DEPENDENCY ${SPHINX_BUILD}/conf.py

    COMMENT "Generating documentation with Sphinx")

add_custom_target(Sphinx DEPENDS ${SPHINX_SOURCE} ${SPHINX_BUILD}/conf.py ${SPHINX_INDEX_FILE})

add_custom_target(documentation ALL DEPENDS Sphinx)

# For MSVC IDE
set_target_properties(
        Sphinx documentation	

		PROPERTIES

		FOLDER DOC
)

# Install documentation
install(
    DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
    DESTINATION ${CMAKE_INSTALL_DATADIR}
)