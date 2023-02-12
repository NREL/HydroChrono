

add_executable(h5fileinfo_t01 h5fileinfo_t01.cpp)
target_link_libraries(h5fileinfo_t01 HydroChrono)

add_executable(chloadaddedmass_t01 chloadaddedmass_t01.cpp)
target_link_libraries(chloadaddedmass_t01 HydroChrono)

add_executable(chrono_error_t01 chrono_error_t01.cpp)
target_link_libraries(chrono_error_t01 HydroChrono)

# ============
# TESTS
# ============
if(TARGET h5fileinfo_t01)
        add_test (
                NAME h5fileinfo_01
                COMMAND $<TARGET_FILE:h5fileinfo_t01> ${HYDROCHRONO_DATA_DIR}
        )
        set_tests_properties(
            h5fileinfo_01
                PROPERTIES LABELS "examples;small;core"
        )
endif(TARGET h5fileinfo_t01)

if(TARGET chloadaddedmass_t01)
        add_test (
                NAME chloadaddedmass_01
                COMMAND $<TARGET_FILE:chloadaddedmass_t01> ${HYDROCHRONO_DATA_DIR}
        )
        set_tests_properties(
                chloadaddedmass_01
                PROPERTIES LABELS "examples;small;core"
        )
endif(TARGET chloadaddedmass_t01)


if(TARGET chrono_error_t01)
        add_test (
                NAME chrono_error_01
                COMMAND $<TARGET_FILE:chrono_error_t01> ${HYDROCHRONO_DATA_DIR}
        )
        set_tests_properties(
                chrono_error_01
                PROPERTIES LABELS "examples;small;core"
        )
endif(TARGET chrono_error_t01)