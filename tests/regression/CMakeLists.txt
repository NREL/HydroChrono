# Sphere Decay Regression Test
# ===========================

# Create executable in model-based directory
add_executable(sphere_decay_test)

set(SPHERE_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sphere)
set(SPHERE_RELEASE_DIR ${CMAKE_BINARY_DIR}/tests/regression/Release/sphere)
set(SPHERE_RESULTS_DIR ${SPHERE_RELEASE_DIR}/results)

# Use model-based source and let simulation create results automatically
set_target_properties(sphere_decay_test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${SPHERE_RELEASE_DIR}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${SPHERE_RELEASE_DIR}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${SPHERE_RELEASE_DIR}
)

target_sources(
    sphere_decay_test
    PRIVATE
        ${SPHERE_SRC_DIR}/demo_sphere_decay.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
    target_compile_definitions(sphere_decay_test 
        PRIVATE 
        HYDROCHRONO_HAVE_IRRLICHT=1
        "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    )
endif()

target_include_directories(
    sphere_decay_test
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../../src
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../../src/gui
)

target_link_libraries(
    sphere_decay_test
    PRIVATE
        HydroChrono
        HydroChronoGUI
)

# Register as CTest test
if(TARGET sphere_decay_test)
    add_test(
        NAME sphere_decay_regression
        COMMAND ${SPHERE_RELEASE_DIR}/sphere_decay_test.exe ${HYDROCHRONO_DATA_DIR} --nogui
    )
    
    set_tests_properties(
        sphere_decay_regression
        PROPERTIES 
        LABELS "regression;sphere;decay;small;core"
        FIXTURES_SETUP sphere_decay_regression_file
        WORKING_DIRECTORY ${SPHERE_RELEASE_DIR}
        ENVIRONMENT "${TEST_ENVIRONMENT}"
    )

    # Reference comparison test - use the new centralized reference data location
    set(FILE_RST ${SPHERE_RESULTS_DIR}/sphere_decay.txt)
    set(FILE_REF ${CMAKE_SOURCE_DIR}/tests/regression/reference_data/sphere/decay/sphere_decay_hc_data.txt)

    add_test(
        NAME sphere_decay_regression_ref
        COMMAND ${Python3_EXECUTABLE} ${SPHERE_SRC_DIR}/compare.py ${FILE_REF} ${FILE_RST}
    )
    
    set_tests_properties(
        sphere_decay_regression_ref
        PROPERTIES 
        LABELS "regression;sphere;decay;reference;core"
        FIXTURES_REQUIRED sphere_decay_regression_file
        WORKING_DIRECTORY ${SPHERE_RESULTS_DIR}
        ENVIRONMENT "${TEST_ENVIRONMENT}"
    )
endif() 