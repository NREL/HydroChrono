# =====================
# SPHERE_DECAY
# =====================
add_executable(demo_sphere_decay)

target_sources(
    demo_sphere_decay

    PRIVATE
        sphere/demo_sphere_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_sphere_decay 
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_sphere_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(demo_sphere_decay 
	PRIVATE
	HydroChronoGUI
	HydroChrono
	
)


# =====================
# SPHERE_REG_WAVES
# =====================
add_executable(demo_sphere_reg_waves)

target_sources(
    demo_sphere_reg_waves

    PRIVATE
        sphere/demo_sphere_reg_waves.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_sphere_reg_waves 
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_sphere_reg_waves

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(demo_sphere_reg_waves 
    PRIVATE 
	HydroChronoGUI	
    HydroChrono
)


# =====================
# SPHERE_IRREG_WAVES
# =====================
add_executable(demo_sphere_irreg_waves)

target_sources(
    demo_sphere_irreg_waves

    PRIVATE
        sphere/demo_sphere_irreg_waves.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_sphere_irreg_waves 
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_sphere_irreg_waves

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(demo_sphere_irreg_waves 
    PRIVATE 
    HydroChronoGUI
    HydroChrono
)

# =====================
# RM3_DECAY
# =====================
add_executable(demo_rm3_decay)

target_sources(
    demo_rm3_decay

    PRIVATE
    rm3/demo_rm3_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_rm3_decay
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_rm3_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(demo_rm3_decay
    PRIVATE 
	HydroChronoGUI	
    HydroChrono
)


# =====================
# RM3_REG_WAVES
# =====================
add_executable(demo_rm3_reg_waves)

target_sources(
    demo_rm3_reg_waves

    PRIVATE
    rm3/demo_rm3_reg_waves.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_rm3_reg_waves
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_rm3_reg_waves

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(demo_rm3_reg_waves 
    PRIVATE
	HydroChronoGUI	
    HydroChrono
)


# =====================
# OSWEC_DECAY
# =====================
add_executable(demo_oswec_decay)

target_sources(
    demo_oswec_decay

    PRIVATE
    oswec/demo_oswec_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_oswec_decay
		PRIVATE 
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_oswec_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(demo_oswec_decay
    PRIVATE
	HydroChronoGUI
    HydroChrono
)


# =====================
# F3OF
# =====================
add_executable(demo_f3of_decay)

target_sources(
    demo_f3of_decay

    PRIVATE
    f3of/demo_F3OF_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_f3of_decay
		PRIVATE
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_f3of_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(demo_f3of_decay
    PRIVATE
	HydroChronoGUI	
    HydroChrono
)


# =====================
# DeepCWind
# =====================
add_executable(demo_DeepCWind_decay)

target_sources(
    demo_DeepCWind_decay

    PRIVATE
    DeepCWind/demo_DeepCWind_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
	target_compile_definitions(demo_DeepCWind_decay
		PRIVATE
		HYDROCHRONO_HAVE_IRRLICHT=1
	)
endif()

target_include_directories(
	demo_DeepCWind_decay

	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(demo_DeepCWind_decay
    PRIVATE
    HydroChrono
)