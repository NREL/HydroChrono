# Manual/Diagnostic Tests for HydroChrono
# =========================================
#
# These are visual verification tools for comparing convolution vs state-space
# radiation damping. They output data files that can be plotted using the
# accompanying Python script.
#
# Usage:
#   cmake --build build --config Release --target sphere_decay_ss_comparison
#   cd build/bin/tests/manual && ./sphere_decay_ss_comparison
#   python plot_model_comparison.py sphere --save
#
#   cmake --build build --config Release --target oswec_decay_ss_comparison  
#   cd build/bin/tests/manual && ./oswec_decay_ss_comparison
#   python plot_model_comparison.py oswec --save
#
# Output structure:
#   build/bin/tests/manual/
#     sphere_decay_ss_comparison.exe
#     oswec_decay_ss_comparison.exe
#     plot_model_comparison.py
#     output/
#       sphere_ss_verification.png
#       oswec_ss_verification.png
#       *_kernel_fit.csv, *_performance.csv, *_convolution.csv, *_state_space.csv

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
# Sphere Decay Comparison: Convolution vs State-Space
# =============================================================================
add_executable(sphere_decay_ss_comparison
    sphere_decay_ss_comparison.cpp
)

target_include_directories(sphere_decay_ss_comparison 
    PRIVATE 
    ${CHRONO_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(sphere_decay_ss_comparison 
    PRIVATE 
    HydroChrono
    ${CHRONO_LIBS}
    Eigen3::Eigen
)

# =============================================================================
# OSWEC Decay Comparison: Convolution vs State-Space (2-body system)
# =============================================================================
add_executable(oswec_decay_ss_comparison
    oswec_decay_ss_comparison.cpp
)

target_include_directories(oswec_decay_ss_comparison 
    PRIVATE 
    ${CHRONO_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(oswec_decay_ss_comparison 
    PRIVATE 
    HydroChrono
    ${CHRONO_LIBS}
    Eigen3::Eigen
)

# Copy the unified plotting script to the build output directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/plot_model_comparison.py
    ${MANUAL_TEST_OUTPUT_DIR}/plot_model_comparison.py
    COPYONLY
)
