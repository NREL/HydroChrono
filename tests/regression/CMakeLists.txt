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
        COMMAND python compare.py ${FILE_REF} ${FILE_RST}
    )
    
    set_tests_properties(
        sphere_decay_regression_ref
        PROPERTIES 
        LABELS "regression;sphere;decay;reference;core"
        FIXTURES_REQUIRED sphere_decay_regression_file
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/regression/sphere
        ENVIRONMENT "PATH=${Python3_ROOT_DIR};$ENV{PATH};HYDROCHRONO_BUILD_DIR=${CMAKE_BINARY_DIR}"
    )
endif() 

# F3OF DT3 Regression Test
# ========================

# Create executable in model-based directory
add_executable(f3of_dt3_test)

set(F3OF_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/f3of)
set(F3OF_RELEASE_DIR ${CMAKE_BINARY_DIR}/tests/regression/Release/f3of)
set(F3OF_RESULTS_DIR ${F3OF_RELEASE_DIR}/results)

# Use model-based source and let simulation create results automatically
set_target_properties(f3of_dt3_test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${F3OF_RELEASE_DIR}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${F3OF_RELEASE_DIR}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${F3OF_RELEASE_DIR}
)

target_sources(
    f3of_dt3_test
    PRIVATE
        ${F3OF_SRC_DIR}/f3of_dt3_test.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
    target_compile_definitions(f3of_dt3_test 
        PRIVATE 
        HYDROCHRONO_HAVE_IRRLICHT=1
        "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    )
endif()

target_include_directories(
    f3of_dt3_test
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../../src
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../../src/gui
)

target_link_libraries(
    f3of_dt3_test
    PRIVATE
        HydroChrono
        HydroChronoGUI
)

# Register as CTest test
if(TARGET f3of_dt3_test)
    # Copy test data to the output directory
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR}/tests/regression/Release)
    
    add_test(
        NAME f3of_dt3_regression
        COMMAND ${F3OF_RELEASE_DIR}/f3of_dt3_test.exe ${HYDROCHRONO_DATA_DIR} --nogui
    )
    
    set_tests_properties(
        f3of_dt3_regression
        PROPERTIES 
        LABELS "regression;f3of;dt3;decay;core"
        FIXTURES_SETUP f3of_dt3_regression_file
        WORKING_DIRECTORY ${F3OF_RELEASE_DIR}
        ENVIRONMENT "${TEST_ENVIRONMENT}"
    )

    # Reference comparison test - use the new centralized reference data location
    set(FILE_RST ${F3OF_RESULTS_DIR}/f3of_dt3.txt)
    set(FILE_REF ${CMAKE_SOURCE_DIR}/tests/regression/reference_data/f3of/dt3/f3of_dt3_hc_data.txt)

    add_test(
        NAME f3of_dt3_regression_ref
        COMMAND python compare.py ${FILE_REF} ${FILE_RST}
    )
    
    set_tests_properties(
        f3of_dt3_regression_ref
        PROPERTIES 
        LABELS "regression;f3of;dt3;decay;reference;core"
        FIXTURES_REQUIRED f3of_dt3_regression_file
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/regression/f3of
        ENVIRONMENT "PATH=${Python3_ROOT_DIR};$ENV{PATH};HYDROCHRONO_BUILD_DIR=${CMAKE_BINARY_DIR}"
    )
endif() 