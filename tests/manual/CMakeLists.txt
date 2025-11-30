# Manual/Diagnostic Tests for HydroChrono
# =========================================
#
# These are NOT automated tests - they are diagnostic tools for visual
# verification and debugging. They output data files that can be plotted
# using the accompanying Python scripts.
#
# Build: cmake --build build --target radiation_ss_visual_check
# Run:   ./radiation_ss_visual_check
# Plot:  python plot_radiation_ss.py
#
# Output structure:
#   build/bin/tests/manual/
#     radiation_ss_visual_check.exe
#     plot_radiation_ss.py
#     output/
#       radiation_ss_*.csv
#       plot_*.png

# Output to tests/manual/ subdirectory
set(MANUAL_TEST_OUTPUT_DIR ${HYDROCHRONO_TEST_OUTPUT_DIR}/manual)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${MANUAL_TEST_OUTPUT_DIR})
if(CMAKE_CONFIGURATION_TYPES)
    foreach(cfg ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER "${cfg}" cfg_uc)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${cfg_uc} ${MANUAL_TEST_OUTPUT_DIR})
    endforeach()
endif()

# =============================================================================
# State-Space Radiation Model Visual Verification
# =============================================================================
add_executable(radiation_ss_visual_check
    radiation_ss_visual_check.cpp
    ../../src/hydro/radiation/radiation_ss_model.cpp
)

target_include_directories(radiation_ss_visual_check 
    PRIVATE 
    ../../src
    ../../include
)

target_link_libraries(radiation_ss_visual_check 
    PRIVATE 
    Eigen3::Eigen
)

# Copy the Python plotting script to the build output directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/plot_radiation_ss.py
    ${MANUAL_TEST_OUTPUT_DIR}/plot_radiation_ss.py
    COPYONLY
)

