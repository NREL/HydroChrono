find_package(Sphinx REQUIRED)

if(NOT Sphinx_FOUND)
    message(WARNING "Can't build programmer's docs. Sphinx not found!")
    message(WARNING "Disable user doc option or provide a sphinx environement")    
    return()
endif(NOT Sphinx_FOUND)

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