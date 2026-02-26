find_package(Python3 COMPONENTS Interpreter REQUIRED)

# ---------------------------------------------------------------------------
# build_verification_test(TEST_GROUP TEST_NAME)
#
# Registers a verification test (HydroChrono vs external reference code).
# Same fixture pattern as build_regression_test() in tests/regression/ but
# uses "verification" CTest labels and reads reference data from
# data/verification/<TEST_GROUP>/.
# ---------------------------------------------------------------------------
function(build_verification_test TEST_GROUP TEST_NAME)

   message(STATUS "...${TEST_NAME} (verification)")

   set(TEST_SRC  ${TEST_GROUP}/test_${TEST_NAME}.cpp)
   set(TEST_EXE  test_${TEST_NAME})
   set(TEST_CMP  ${TEST_NAME}_verification)
   set(COMPARE_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_GROUP}/compare_${TEST_NAME}.py)
   set(FIXTURE_FILE   ${TEST_NAME}_verification_file)

   # Reference files live under data/verification/<TEST_GROUP>/reference/
   set(FILE_REF_DIR ${CMAKE_SOURCE_DIR}/data/verification/${TEST_GROUP}/reference)
   set(FILE_RST     results/tests/${TEST_GROUP}/results_${TEST_NAME}.txt)

   add_executable(${TEST_EXE} ${TEST_SRC})
   source_group("" FILES ${TEST_SRC})

   target_compile_definitions(${TEST_EXE}
       PRIVATE
       "REFERENCE_DIR=\"${FILE_REF_DIR}\""
       "RESULTS_DIR_NAME=\"${TEST_GROUP}\""
       "RESULTS_FILE_NAME=\"results_${TEST_NAME}\""
   )

   if(MSVC)
       set_target_properties(${TEST_EXE} PROPERTIES
           VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${TEST_EXE}>")
   endif()

   target_link_libraries(${TEST_EXE} PRIVATE HydroChrono HydroChronoGUI)

   # ── Execute test ──────────────────────────────────────────────────────
   add_test(
       NAME ${TEST_EXE}
       WORKING_DIRECTORY "$<TARGET_FILE_DIR:${TEST_EXE}>"
       COMMAND $<TARGET_FILE:${TEST_EXE}> --no_gui --no_plot)

   set_tests_properties(
       ${TEST_EXE}
       PROPERTIES
       LABELS "verification;${TEST_GROUP};execute"
       FIXTURES_SETUP ${FIXTURE_FILE}
       ENVIRONMENT "HYDROCHRONO_DATA_DIR=${HC_BUILD_DATA}"
   )

   # ── Comparison test ───────────────────────────────────────────────────
   add_test(
       NAME ${TEST_CMP}
       WORKING_DIRECTORY "$<TARGET_FILE_DIR:${TEST_EXE}>"
       COMMAND ${Python3_EXECUTABLE} ${COMPARE_SCRIPT}
               ${FILE_REF_DIR}
               $<TARGET_FILE_DIR:${TEST_EXE}>/${FILE_RST}
   )

   if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.22")
       set_tests_properties(
           ${TEST_CMP}
           PROPERTIES
           LABELS "verification;${TEST_GROUP};reference"
           FIXTURES_REQUIRED ${FIXTURE_FILE}
           ENVIRONMENT_MODIFICATION
               "PATH=path_list_prepend:${Python3_ROOT_DIR};PATH=path_list_prepend:${Python3_ROOT_DIR}/Library/bin;PATH=path_list_prepend:${Python3_ROOT_DIR}/DLLs"
           ENVIRONMENT "HYDROCHRONO_BUILD_DIR=${CMAKE_BINARY_DIR}"
       )
   else()
       set_tests_properties(
           ${TEST_CMP}
           PROPERTIES
           LABELS "verification;${TEST_GROUP};reference"
           FIXTURES_REQUIRED ${FIXTURE_FILE}
           ENVIRONMENT "PATH=${Python3_ROOT_DIR};HYDROCHRONO_BUILD_DIR=${CMAKE_BINARY_DIR}"
       )
   endif()

endfunction()

# ===========================================================================
# Register verification tests
# ===========================================================================

message(STATUS "Add HydroChrono verification tests")

if(HYDROCHRONO_ENABLE_MOORDYN)
    build_verification_test(rm3_mooring  rm3_mooring)
else()
    message(STATUS "  (none -- rm3_mooring requires HYDROCHRONO_ENABLE_MOORDYN)")
endif()
