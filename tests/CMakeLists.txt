

add_executable(h5fileinfo_t01 h5fileinfo_t01.cpp)


target_link_libraries(h5fileinfo_t01 HydroChrono)


if(TARGET h5fileinfo_t01)
        add_test (
                NAME h5fileinfo_01
                COMMAND $<TARGET_FILE:h5fileinfo_t01> ${HYDRO_CHRONO_DATA_DIR}
        )
        set_tests_properties(
            h5fileinfo_01
                PROPERTIES LABELS "examples;small;core"
        )
endif(TARGET h5fileinfo_t01)
