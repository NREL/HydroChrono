# =====================
# SPHERE_DECAY
# =====================
add_executable(sphere_decay)

target_sources(
    sphere_decay

    PRIVATE
        sphere/sphere_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(sphere_decay 
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	sphere_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(sphere_decay 
	PRIVATE
	HydroChronoGUI
	HydroChrono
	
)


# =====================
# SPHERE_REG_WAVES
# =====================
add_executable(sphere_reg_waves)

target_sources(
    sphere_reg_waves

    PRIVATE
        sphere/sphere_reg_waves.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(sphere_reg_waves 
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	sphere_reg_waves

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(sphere_reg_waves 
    PRIVATE 
	HydroChronoGUI	
    HydroChrono
)


# =====================
# RM3_DECAY
# =====================
add_executable(rm3_decay)

target_sources(
    rm3_decay

    PRIVATE
    rm3/rm3_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(rm3_decay
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	rm3_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(rm3_decay
    PRIVATE 
	HydroChronoGUI	
    HydroChrono
)


# =====================
# RM3_REG_WAVES
# =====================
add_executable(rm3_reg_waves)

target_sources(
    rm3_reg_waves

    PRIVATE
    rm3/rm3_reg_waves.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(rm3_reg_waves
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	rm3_reg_waves

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(rm3_reg_waves 
    PRIVATE
	HydroChronoGUI	
    HydroChrono
)


# =====================
# OSWEC_DECAY
# =====================
add_executable(oswec_decay)

target_sources(
    oswec_decay

    PRIVATE
    oswec/oswec_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(oswec_decay
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	oswec_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(oswec_decay
    PRIVATE
    HydroChrono
)


# =====================
# F3OF
# =====================
add_executable(f3of)

target_sources(
    f3of

    PRIVATE
    f3of/F3OF.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(f3of
		PRIVATE
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	f3of

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(f3of
    PRIVATE
    HydroChrono
)


# =====================
# DeepCWind
# =====================
add_executable(demo_DeepCWind)

target_sources(
    demo_DeepCWind

    PRIVATE
    DeepCWind/demo_DeepCWind.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_DeepCWind
		PRIVATE
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_DeepCWind

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(demo_DeepCWind
    PRIVATE
    HydroChrono
)