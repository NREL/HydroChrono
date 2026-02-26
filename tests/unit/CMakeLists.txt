# ==============================================================================
# Unit tests — standalone executables, no external data dependencies
# ==============================================================================

# Helper: build a self-contained unit test and register with CTest
function(build_unit_test TEST_NAME)
    set(TEST_SRC ${TEST_NAME}.cpp)
    set(TEST_EXE ${TEST_NAME})

    add_executable(${TEST_EXE} ${TEST_SRC})
    target_link_libraries(${TEST_EXE} PRIVATE HydroChrono)
    target_compile_features(${TEST_EXE} PRIVATE ${HC_CXX_STANDARD})

    add_test(
        NAME ${TEST_EXE}
        WORKING_DIRECTORY "$<TARGET_FILE_DIR:${TEST_EXE}>"
        COMMAND $<TARGET_FILE:${TEST_EXE}>
    )

    set_tests_properties(${TEST_EXE}
        PROPERTIES LABELS "unit"
    )

    # DLL search paths on Windows (reuse parent-scope variable if set)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows" AND DEFINED TEST_ENVIRONMENT)
        set_tests_properties(${TEST_EXE}
            PROPERTIES ENVIRONMENT "${TEST_ENVIRONMENT}"
        )
    endif()
endfunction()

# --- Register unit tests ---
message(STATUS "Add HydroChrono unit tests")
if(CHRONO_HAS_LOAD_HYDRODYNAMICS)
    build_unit_test(test_added_mass_determinism)
else()
    message(STATUS "  Skipping test_added_mass_determinism (requires ChLoadHydrodynamics)")
endif()
build_unit_test(test_radiation_ss_model)
build_unit_test(test_radiation_ss_fitter)
