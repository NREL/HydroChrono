

add_executable(sphere_decay sphere/sphere_decay.cpp)
if(HYDRO_CHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(sphere_decay 
		PRIVATE 
		HYDRO_CHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	sphere_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(sphere_decay HydroChrono)


# add_executable(sphere_reg_waves "demos/sphere/sphere_reg_waves.cpp")
# add_executable(rm3_decay "demos/rm3/rm3_decay.cpp")
# add_executable(rm3_reg_waves "demos/rm3/rm3_reg_waves.cpp")
# add_executable(oswec_decay "demos/oswec/oswec_decay.cpp")
#SWIG_ADD_LIBRARY(py_wec_chrono LANGUAGE python SOURCES hydro_forces_no_viz.i)



#--------------------------------------------------------------
# Set properties for your executable target
# 
# Note that here we define a macro CHRONO_DATA_DIR which will
# contain the path to the Chrono data directory, either in its
# source tree (if using a build version of Chrono), or in its
# install tree (if using an installed version of Chrono).
#--------------------------------------------------------------



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



# target_link_libraries(sphere_reg_waves HydroChrono)
# target_link_libraries(rm3_decay HydroChrono)
# target_link_libraries(rm3_reg_waves HydroChrono)
# target_link_libraries(oswec_decay HydroChrono)
#SWIG_LINK_LIBRARIES(py_wec_chrono wec_chrono ${PYTHON_LIBRARIES})
