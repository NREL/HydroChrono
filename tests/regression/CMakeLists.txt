find_package(Python3 COMPONENTS Interpreter REQUIRED)

function(build_regression_test TEST_GROUP TEST_NAME MULTIPLE)

   message(STATUS "...${TEST_NAME}")
   
   set(TEST_SRC ${TEST_GROUP}/test_${TEST_NAME}.cpp)
   set(TEST_EXE test_${TEST_NAME})
   set(TEST_REF ${TEST_NAME}_regression)
   set(COMPARE_SCRIPT ${CMAKE_SOURCE_DIR}/tests/regression/${TEST_GROUP}/compare_${TEST_NAME}.py)
   set(FIXTURE_FILE ${TEST_NAME}_regression_file)
   if(MULTIPLE)
     set(FILE_REF ${CMAKE_SOURCE_DIR}/data/reference_data/${TEST_GROUP}/hc_ref_${TEST_NAME}_1.txt)
     set(FILE_RST results/tests/${TEST_GROUP}/results_${TEST_NAME}_1.txt)
   else()
     set(FILE_REF ${CMAKE_SOURCE_DIR}/data/reference_data/${TEST_GROUP}/hc_ref_${TEST_NAME}.txt)
     set(FILE_RST results/tests/${TEST_GROUP}/results_${TEST_NAME}.txt)
   endif()

   add_executable(${TEST_EXE} ${TEST_SRC})
   source_group("" FILES ${TEST_SRC})
      
   target_compile_definitions(${TEST_EXE} 
	   PRIVATE
       "REFERENCE_FILE_NAME=\"${FILE_REF}\""
       "RESULTS_DIR_NAME=\"${TEST_GROUP}\""
       "RESULTS_FILE_NAME=\"results_${TEST_NAME}\""
   )

   if(MSVC)
       set_target_properties(${TEST_EXE} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${TEST_EXE}>")
   endif()

   target_link_libraries(${TEST_EXE} PRIVATE HydroChrono HydroChronoGUI)

   add_test(
       NAME ${TEST_EXE}
       WORKING_DIRECTORY "$<TARGET_FILE_DIR:${TEST_EXE}>"
       COMMAND $<TARGET_FILE:${TEST_EXE}> --no_gui --no_plot)

   # Inject HYDROCHRONO_DATA_DIR so tests find demo assets (meshes, H5 files, etc.)
   # regardless of what the user's shell environment has set.
   set_tests_properties(
       ${TEST_EXE}
       PROPERTIES 
       LABELS "regression;${TEST_GROUP};execute"
       FIXTURES_SETUP ${FIXTURE_FILE}
       ENVIRONMENT "HYDROCHRONO_DATA_DIR=${HC_BUILD_DATA}"
   )

   # Reference comparison test - use the new centralized reference data location

   add_test(
       NAME ${TEST_REF}
       WORKING_DIRECTORY "$<TARGET_FILE_DIR:${TEST_EXE}>"
       COMMAND ${Python3_EXECUTABLE} ${COMPARE_SCRIPT} ${FILE_REF} $<TARGET_FILE_DIR:${TEST_EXE}>/${FILE_RST}
   )
   
   if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.22")
       set_tests_properties(
           ${TEST_REF}
           PROPERTIES
           LABELS "regression;${TEST_GROUP};reference"
           FIXTURES_REQUIRED ${FIXTURE_FILE}
           ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${Python3_ROOT_DIR};PATH=path_list_prepend:${Python3_ROOT_DIR}/Library/bin;PATH=path_list_prepend:${Python3_ROOT_DIR}/DLLs"
           ENVIRONMENT "HYDROCHRONO_BUILD_DIR=${CMAKE_BINARY_DIR}"
       )
   else()
       # CMake < 3.22 lacks ENVIRONMENT_MODIFICATION; set only the Python root in
       # PATH (semicolons inside values can't be escaped).  Users on older CMake
       # should activate their Python environment before running ctest.
       set_tests_properties(
           ${TEST_REF}
           PROPERTIES
           LABELS "regression;${TEST_GROUP};reference"
           FIXTURES_REQUIRED ${FIXTURE_FILE}
           ENVIRONMENT "PATH=${Python3_ROOT_DIR};HYDROCHRONO_BUILD_DIR=${CMAKE_BINARY_DIR}"
       )
   endif()

endfunction()

# Build tests
# ===========

message(STATUS "Add HydroChrono regression tests")

build_regression_test(sphere  sphere_decay            FALSE)
build_regression_test(sphere  sphere_decay_ss         FALSE)
build_regression_test(sphere  sphere_reg_waves        TRUE)
build_regression_test(sphere  sphere_irreg_waves      FALSE)
build_regression_test(sphere  sphere_irreg_waves_ss   FALSE)
build_regression_test(sphere  sphere_irreg_waves_eta  FALSE)
build_regression_test(sphere  sphere_irreg_waves_eta_consistency  FALSE)
build_regression_test(f3of    f3of_dt1                FALSE)
build_regression_test(f3of    f3of_dt2                FALSE)
build_regression_test(f3of    f3of_dt3                FALSE)
build_regression_test(oswec   oswec_decay             FALSE)
build_regression_test(oswec   oswec_decay_ss          FALSE)
build_regression_test(oswec   oswec_irreg_waves       FALSE)
build_regression_test(oswec   oswec_irreg_waves_ss    FALSE)
build_regression_test(oswec   oswec_reg_waves         TRUE)
build_regression_test(rm3     rm3_decay               FALSE)
build_regression_test(rm3     rm3_reg_waves           FALSE)

# Regression Test Report Generation
# =================================

# Create documentation directory next to other regression artefacts
set(REGRESSION_DOCS_DIR ${CMAKE_BINARY_DIR}/bin/report)
file(MAKE_DIRECTORY ${REGRESSION_DOCS_DIR})

# Add a test that generates the regression report (markdown + PDF) after all tests complete
add_test(
    NAME regression_report_generation
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/regression/utilities/generate_report.py --build-dir ${CMAKE_BINARY_DIR} --output-dir ${REGRESSION_DOCS_DIR} --pdf
)

# Make sure report generation runs after all other regression tests
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.22")
    set_tests_properties(
        regression_report_generation
        PROPERTIES
        LABELS "regression;report;documentation;core"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${Python3_ROOT_DIR};PATH=path_list_prepend:${Python3_ROOT_DIR}/Library/bin;PATH=path_list_prepend:${Python3_ROOT_DIR}/DLLs"
        ENVIRONMENT "HYDROCHRONO_BUILD_DIR=${CMAKE_BINARY_DIR}"
        DEPENDS "sphere_decay_regression;sphere_decay_ss_regression;sphere_reg_waves_regression;sphere_irreg_waves_regression;sphere_irreg_waves_ss_regression;sphere_irreg_waves_eta_regression;sphere_irreg_waves_eta_consistency_regression;f3of_dt1_regression;f3of_dt2_regression;f3of_dt3_regression;oswec_decay_regression;oswec_decay_ss_regression;oswec_irreg_waves_regression;oswec_irreg_waves_ss_regression;oswec_reg_waves_regression;rm3_decay_regression;rm3_reg_waves_regression"
    )
else()
    set_tests_properties(
        regression_report_generation
        PROPERTIES
        LABELS "regression;report;documentation;core"
        ENVIRONMENT "PATH=${Python3_ROOT_DIR};HYDROCHRONO_BUILD_DIR=${CMAKE_BINARY_DIR}"
        DEPENDS "sphere_decay_regression;sphere_decay_ss_regression;sphere_reg_waves_regression;sphere_irreg_waves_regression;sphere_irreg_waves_ss_regression;sphere_irreg_waves_eta_regression;sphere_irreg_waves_eta_consistency_regression;f3of_dt1_regression;f3of_dt2_regression;f3of_dt3_regression;oswec_decay_regression;oswec_decay_ss_regression;oswec_irreg_waves_regression;oswec_irreg_waves_ss_regression;oswec_reg_waves_regression;rm3_decay_regression;rm3_reg_waves_regression"
    )
endif()
