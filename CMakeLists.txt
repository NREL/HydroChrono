# ===============================================================================
# ------ Prerequisites & Policies -----------------------------------------------
# ===============================================================================

# cmake_minimum_required(3.18.2): need CMP0091 policy + target properties support
cmake_minimum_required(VERSION 3.18.2)

# CMP0091: Enable MSVC runtime library selection via CMAKE_MSVC_RUNTIME_LIBRARY
cmake_policy(SET CMP0091 NEW)
# CMP0074: Use <PackageName>_ROOT variables for find_package()
cmake_policy(SET CMP0074 NEW)

# Guard against in-source builds to prevent polluting source directory
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif()

# ===============================================================================
# ------- Project Definition & Defaults -----------------------------------------
# ===============================================================================

project(HydroChrono 
    VERSION 0.5.0
    DESCRIPTION "Hydrodynamics for Project Chrono."
    LANGUAGES CXX
)

# MSVC: CMake may pass "SYSTEM" include directories as external headers (/external:I).
# Some third-party headers rely on nested quoted includes that can fail intermittently when treated as external headers.
# Force SYSTEM include dirs to be emitted as normal /I include paths instead.
if(MSVC)
    set(CMAKE_INCLUDE_SYSTEM_FLAG_C "/I")
    set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "/I")
endif()

# ===============================================================================
# -------- Project Metadata ----------------------------------------------------
# ===============================================================================

# Project metadata derived from project() command - single source of truth
set(HYDROCHRONO_NAME "${PROJECT_NAME}")
set(HYDROCHRONO_VERSION "${PROJECT_VERSION}")
set(HYDROCHRONO_DESCRIPTION "${PROJECT_DESCRIPTION}")

message(STATUS "Project: ${HYDROCHRONO_NAME} v${HYDROCHRONO_VERSION}")
message(STATUS "Description: ${HYDROCHRONO_DESCRIPTION}")

# Generate version.h header from template
configure_file(${CMAKE_SOURCE_DIR}/cmake/version.h.in
               ${CMAKE_BINARY_DIR}/hydroc/version.h @ONLY)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Deterministic dependency discovery (avoid ambient PATH/registry drift)
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
set(CMAKE_FIND_USE_PACKAGE_REGISTRY OFF)
set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH OFF)

# Force Release by default so CI builds are optimized and users get good performance
if(NOT DEFINED HYDROCHRONO_DEFAULT_BUILD_TYPE)
    set(HYDROCHRONO_DEFAULT_BUILD_TYPE "Release")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${HYDROCHRONO_DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE ${HYDROCHRONO_DEFAULT_BUILD_TYPE}
        CACHE STRING "Choose the type of build." FORCE)
    mark_as_advanced(CMAKE_BUILD_TYPE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ===============================================================================
# ------ User Options -----------------------------------------------------------
# ===============================================================================

option(HYDROCHRONO_ENABLE_TESTS "Enable tests" ON)
option(HYDROCHRONO_ENABLE_IRRLICHT "Enable Irrlicht visualization library" OFF)
option(HYDROCHRONO_ENABLE_VSG "Enable VSG visualization library" OFF)
option(HYDROCHRONO_ENABLE_DEMOS "Enable demo executables" OFF)
option(HYDROCHRONO_ENABLE_YAML_RUNNER "Enable YAML-based CLI runner" OFF)

# ===============================================================================
# ------ Build setup ------------------------------------------------------------
# ===============================================================================

message(STATUS "")

# Fix for VS 2017 15.8 and newer to handle alignment specification with Eigen
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if(MSVC AND ${MSVC_VERSION} GREATER_EQUAL 1915)
        add_compile_definitions("ENABLE_EXTENDED_ALIGNED_STORAGE")
    endif()
endif()

# Disable common MSVC warnings that pollute build output
if(MSVC)
    add_definitions(
        -D_CRT_SECURE_NO_DEPRECATE   # Disable CRT deprecation warnings
        -D_SCL_SECURE_NO_DEPRECATE   # Disable SCL deprecation warnings
        -DNOMINMAX                   # Prevent Windows.h from defining min/max macros
    )
    add_compile_options(/wd4275)     # Disable non-DLL class used as base for DLL-interface class
    add_compile_options(/wd4251)     # Disable "class needs to have dll-interface" warnings
endif()

# -- HydroChrono data directory

message(STATUS "HydroChrono data directory")

# Destination directory for HydroChrono data in BUILD tree
set(HC_BUILD_DATA "${CMAKE_BINARY_DIR}/data")

# Destination directory for HydroChrono data in INSTALL tree
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(HC_INSTALL_DATA "data")
    set(HC_INSTALL_DEMO "bin")
    set(HC_INSTALL_PYTHON "bin")
else()
    set(HC_INSTALL_DATA "share/chrono/data")
    set(HC_INSTALL_DEMO "share/chrono/bin")
    set(HC_INSTALL_PYTHON "share/chrono/python")
endif()

# Copy data directory to BUILD tree
file(COPY ${CMAKE_SOURCE_DIR}/data/ DESTINATION ${HC_BUILD_DATA})
message(STATUS "  copied to: ${HC_BUILD_DATA}/")

# Install data directory (trailing slash copies contents, not the directory itself)
install(DIRECTORY ${HC_BUILD_DATA}/ DESTINATION ${HC_INSTALL_DATA})
message(STATUS "  installed to: ${CMAKE_INSTALL_PREFIX}/${HC_INSTALL_DATA}/")

# -- Output Directory Structure --

# Centralize all binaries for easier deployment and testing
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Replicate for multi-config generators (Visual Studio, Xcode, Ninja Multi-Config)
if(CMAKE_CONFIGURATION_TYPES)
    foreach(cfg ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER "${cfg}" cfg_uc)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${cfg_uc} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${cfg})
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${cfg_uc} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${cfg})
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${cfg_uc} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${cfg})
    endforeach()
elseif(CMAKE_BUILD_TYPE)
    # Single-config generator with explicit build type
    string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_uc)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${build_type_uc} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${build_type_uc} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${build_type_uc} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE})
endif()

# ===============================================================================
# ------- Find Dependencies -----------------------------------------------------
# ===============================================================================

# -- Chrono package --

message(STATUS "\n---- Find Chrono")

if(NOT Chrono_DIR)
    message("Specify the directory containing a CMake configuration file for Chrono")
    set(Chrono_DIR "" CACHE PATH "The directory containing a CMake configuration file for Chrono.")
    return()
endif()

find_package(Chrono
             CONFIG REQUIRED
             COMPONENTS Parsers)

message(STATUS "----\n")

message(STATUS "Chrono targets: ${CHRONO_TARGETS}")

# Verify required Chrono targets are available
if(NOT TARGET Chrono::Chrono_core)
    message(FATAL_ERROR "Chrono::Chrono_core target not found. Ensure Chrono was built with modern CMake target exports.")
endif()
if(NOT TARGET Chrono::Chrono_parsers)
    message(FATAL_ERROR "Chrono::Chrono_parsers target not found. Ensure Chrono was built with modern CMake target exports.")
endif()

# Ensure we have a compatible Chrono version (check after find_package)
if(DEFINED CHRONO_VERSION)
    if(CHRONO_VERSION VERSION_LESS "9.0.1")
        message(WARNING "HydroChrono now requires Chrono ≥ 9.0.1 for modular target support. Found version: ${CHRONO_VERSION}")
    endif()
else()
    # If CHRONO_VERSION is not defined, assume it's a compatible version
    # since the modern targets are available
    message(STATUS "Chrono version not explicitly set, but modern targets are available - assuming compatible version")
endif()

# To ensure ABI compatibility, use the same C++ standard as the one used to build Chrono
if (CHRONO_CXX_STANDARD)
  set(HC_CXX_STANDARD ${CHRONO_CXX_STANDARD})
  message(STATUS "Set C++ standard same as Chrono: ${HC_CXX_STANDARD}")
else()
  set(HC_CXX_STANDARD cxx_std_17)
  message(STATUS "Set C++ standard to: ${HC_CXX_STANDARD}")
endif()

# Enable Irrlicht support if requested and available
if(HYDROCHRONO_ENABLE_IRRLICHT AND NOT TARGET Chrono::Chrono_irrlicht)
  message(FATAL_ERROR "HYDROCHRONO_ENABLE_IRRLICHT is ON but Chrono::Chrono_irrlicht target not found. Ensure Chrono was built with Irrlicht support.")
endif()

# Enable VSG support if requested and available
if(HYDROCHRONO_ENABLE_VSG AND NOT TARGET Chrono::Chrono_vsg)
  message(FATAL_ERROR "HYDROCHRONO_ENABLE_VSG is ON but Chrono::Chrono_vsg target not found. Ensure Chrono was built with VSG support.")
endif()

# If on a Windows platform, copy the Chrono DLLs to the binary directory
add_CHRONO_DLLS_copy_command()

# Copy selected Chrono data directories to the HydroChrono data directory
if(EXISTS "${CHRONO_DATA_DIR}")
    if(EXISTS "${CHRONO_DATA_DIR}/colormaps")
        file(COPY ${CHRONO_DATA_DIR}/colormaps DESTINATION ${HC_BUILD_DATA}/chrono)
    endif()

    if(EXISTS "${CHRONO_DATA_DIR}/skybox")
        file(COPY ${CHRONO_DATA_DIR}/skybox DESTINATION ${HC_BUILD_DATA}/chrono)
    endif()

    if(EXISTS "${CHRONO_DATA_DIR}/vsg")
        file(COPY ${CHRONO_DATA_DIR}/vsg DESTINATION ${HC_BUILD_DATA}/chrono)
    endif()

    if(EXISTS "${CHRONO_DATA_DIR}/logo_chrono_alpha.png")
        file(COPY ${CHRONO_DATA_DIR}/logo_chrono_alpha.png DESTINATION ${HC_BUILD_DATA}/chrono)
    endif()
endif()

# -- OpenMP support --
find_package(OpenMP REQUIRED COMPONENTS CXX)

# -- HDF5 Integration --
find_package(HDF5 REQUIRED COMPONENTS CXX)

# ===============================================================================
# ----- Create HydroChrono configuration header ---------------------------------
# ===============================================================================

message(STATUS "Generate HydroChrono configuration headers")

set(CH_DATA_DIR "#define CHRONO_DATA_DIR \"${CHRONO_DATA_DIR}\"")

if(HYDROCHRONO_ENABLE_IRRLICHT)
  set(HC_HAS_IRRLICHT "#define HYDROCHRONO_HAVE_IRRLICHT")
else()
  set(HC_HAS_IRRLICHT "#undef HYDROCHRONO_HAVE_IRRLICHT")
endif()

if(HYDROCHRONO_ENABLE_VSG)
  set(HC_HAS_VSG "#define HYDROCHRONO_HAVE_VSG")
else()
  set(HC_HAS_VSG "#undef HYDROCHRONO_HAVE_VSG")
endif()

# For BUILD tree
set(HC_DATA_DIR "#define HC_DATA_DIR \"${HC_BUILD_DATA}\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in
               ${PROJECT_BINARY_DIR}/hydroc/config.h
               @ONLY)

# For INSTALL tree
set(HC_DATA_DIR "#define HC_DATA_DIR \"${CMAKE_INSTALL_PREFIX}/${HC_INSTALL_DATA}\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in
               ${PROJECT_BINARY_DIR}/hydroc/config.h-install
               @ONLY)

# ===============================================================================
# ------- Core Library Target ---------------------------------------------------
# ===============================================================================

list(APPEND HC_INCLUDES_BUILD ${PROJECT_SOURCE_DIR}/include)
list(APPEND HC_INCLUDES_BUILD ${PROJECT_BINARY_DIR})

# -- Main HydroChrono Library --
set(HC_SOURCES
    # IO module (consolidated under src/hydro/io/)
    src/hydro/io/h5_reader.cpp
    src/hydro/io/h5_writer.cpp
    src/hydro/io/simulation_export.cpp
    # Config/YAML module (consolidated under src/hydro/config/)
    src/hydro/config/yaml_parser.cpp
    src/hydro/config/config_loader.cpp
    src/hydro/config/setup_from_yaml.cpp
    src/hydro/config/yaml_discovery.cpp
    # Core facade and utilities
    src/hydro/hydro_system.cpp
    src/hydro/utils/helper.cpp
    # Hydro modules
    src/hydro/waves/wave_base.cpp
    src/hydro/waves/wave_utilities.cpp
    src/hydro/waves/regular_wave.cpp
    src/hydro/waves/irregular_wave.cpp
    src/hydro/radiation/radiation_rirf_processing.cpp
    src/hydro/radiation/radiation_rirf_convolution.cpp
    src/hydro/core/hydro_forces.cpp
    # Chrono coupling layer (bridges Chrono physics with Chrono-free core)
    src/hydro/chrono/added_mass.cpp
    src/hydro/chrono/chrono_state_utils.cpp
    src/hydro/chrono/chrono_hydro_coupler.cpp
    src/hydro/force_components/hydrostatics_component.cpp
    src/hydro/force_components/radiation_component.cpp
    src/hydro/force_components/excitation_component.cpp
    # Logging module (consolidated under src/hydro/logging/)
    src/hydro/logging/logger_backend.cpp
    src/hydro/logging/logging.cpp
)

set(HC_HEADERS
    ${PROJECT_BINARY_DIR}/hydroc/config.h
    ${PROJECT_BINARY_DIR}/hydroc/version.h
    include/hydroc/hydro_system.h
    include/hydroc/helper.h
    include/hydroc/hydro_forces.h
    include/hydroc/logging.h
    include/hydroc/config/hydro_config.h
    include/hydroc/core/force_component.h
    include/hydroc/core/hydro_forces.h
    include/hydroc/core/hydro_types.h
    include/hydroc/core/system_state.h
    include/hydroc/coupling/added_mass.h
    include/hydroc/coupling/chrono_coupler.h
    include/hydroc/io/h5_reader.h
    include/hydroc/io/h5_writer.h
    include/hydroc/io/simulation_export.h
    include/hydroc/radiation/radiation_types.h
    include/hydroc/waves/irregular_wave.h
    include/hydroc/waves/regular_wave.h
    include/hydroc/waves/wave_base.h
)

# Create the library target and configure it
add_library(HydroChrono ${HC_SOURCES} ${HC_HEADERS})

target_compile_features(HydroChrono PUBLIC ${HC_CXX_STANDARD})
target_include_directories(HydroChrono
    PUBLIC
        "$<BUILD_INTERFACE:${HC_INCLUDES_BUILD}>"
        "$<INSTALL_INTERFACE:include>"
)
set_target_properties(HydroChrono PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(HydroChrono PUBLIC Chrono::Chrono_core)

# Add OpenMP
target_link_libraries(HydroChrono PUBLIC OpenMP::OpenMP_CXX)

# Add HDF5 (assume first latest HDF5 version 1.14.6 and prefer static library)
if(TARGET hdf5::hdf5-static)
  message(STATUS "Linking static libraries hdf5::hdf5-static")
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5-static)
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5_cpp-static)
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5_hl-static)
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5_hl_cpp-static)
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5_tools-static)
elseif(TARGET hdf5::hdf5-shared)
  message(STATUS "Linking shared libraries hdf5::hdf5-shared")
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5-shared)
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5_cpp-shared)
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5_hl-shared)
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5_hl_cpp-shared)
  target_link_libraries(HydroChrono PUBLIC hdf5::hdf5_tools-shared)
  list(APPEND HDF5_TARGETS hdf5::hdf5-shared)
  list(APPEND HDF5_TARGETS hdf5::hdf5_cpp-shared)
  list(APPEND HDF5_TARGETS hdf5::hdf5_hl-shared)
  list(APPEND HDF5_TARGETS hdf5::hdf5_hl_cpp-shared)
  list(APPEND HDF5_TARGETS hdf5::hdf5_tools-shared)
elseif(TARGET HDF5::HDF5)
  message(STATUS "Linking HDF5::HDF5")
  target_link_libraries(HydroChrono PUBLIC HDF5::HDF5)
  list(APPEND HDF5_TARGETS HDF5::HDF5)
else()
  message(STATUS "No HDF5 target found")
  target_compile_definitions(HydroChrono PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${HDF5_COMPILE_DEFS}>)
  target_include_directories(HydroChrono PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(HydroChrono PUBLIC ${HDF5_CXX_LIBRARIES})
endif()

target_include_directories(HydroChrono
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Library-specific compile definitions
target_compile_definitions(HydroChrono
    PUBLIC
        "HYDROCHRONO_VERSION=\"${PROJECT_VERSION}\""
        "HYDROCHRONO_BUILD_TYPE=\"$<IF:$<CONFIG:Debug>,Debug,$<IF:$<CONFIG:Release>,Release,$<IF:$<CONFIG:RelWithDebInfo>,RelWithDebInfo,$<IF:$<CONFIG:MinSizeRel>,MinSizeRel,Unknown>>>>\""
    PRIVATE 
        $<$<BOOL:${HYDROCHRONO_ENABLE_LOGGING}>:HYDROCHRONO_ENABLE_LOGGING="1">
)

# ===============================================================================
# ------- Visualization library -------------------------------------------------
# ===============================================================================

set(HCGUI_HEADERS
    include/hydroc/gui/guihelper.h
    src/gui/guihelper_impl.h
)

set(HCGUI_SOURCES
    src/gui/guihelper.cpp
)

if(HYDROCHRONO_ENABLE_IRRLICHT)
  set(HCGUI_SOURCES ${HCGUI_SOURCES}
      src/gui/guihelperIRR.cpp
  )
endif()

if(HYDROCHRONO_ENABLE_VSG)
  set(HCGUI_SOURCES ${HCGUI_SOURCES}
      src/gui/guihelperVSG.cpp
      src/gui/vsg_gui_component.cpp
      src/gui/vsg_lighting.cpp
      src/gui/vsg_materials.cpp
      src/gui/vsg_water_surface.cpp
      src/gui/vsg_radiation_surface.cpp
  )
  set(HCGUI_HEADERS ${HCGUI_HEADERS}
      src/gui/vsg_config.h
      src/gui/vsg_gui_component.h
      src/gui/vsg_lighting.h
      src/gui/vsg_materials.h
      src/gui/vsg_water_surface.h
      src/gui/vsg_radiation_surface.h
  )
endif()

add_library(HydroChronoGUI ${HCGUI_SOURCES} ${HCGUI_HEADERS})

target_compile_features(HydroChronoGUI PUBLIC ${HC_CXX_STANDARD})
target_include_directories(HydroChronoGUI
    PUBLIC
        "$<BUILD_INTERFACE:${HC_INCLUDES_BUILD}>"
        "$<INSTALL_INTERFACE:include>"
)

set_target_properties(HydroChronoGUI PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(HydroChronoGUI
    PUBLIC
        Chrono::Chrono_core
        # Conditionally link Irrlicht and VSG components if enabled
        $<$<BOOL:${HYDROCHRONO_ENABLE_IRRLICHT}>:Chrono::Chrono_irrlicht>
        $<$<BOOL:${HYDROCHRONO_ENABLE_VSG}>:Chrono::Chrono_vsg>
)

# ===============================================================================
# ------- Auxiliary Targets (app, tests, demos) ---------------------------------
# ===============================================================================

# Helper function to configure test environment (Windows DLL paths, etc.)
function(configure_test_environment)
    message(STATUS "Configure test environment")
    
    # HydroChrono DLL paths
    set(HC_DLL_DIRS "${CMAKE_BINARY_DIR}/bin/Release;${CMAKE_BINARY_DIR}/bin/RelWithDebInfo;${CMAKE_BINARY_DIR}/bin/MinSizeRel;${CMAKE_BINARY_DIR}/bin/Debug")
    message(STATUS "  HydroChrono DLL directories:      ${HC_DLL_DIRS}")
    list(APPEND DLL_DIRS ${HC_DLL_DIRS})

    # Chrono DLL paths
    get_target_property(CHRONO_CORE_REL_DLL Chrono::Chrono_core IMPORTED_LOCATION_RELEASE)
    if(CHRONO_CORE_REL_DLL)
        get_filename_component(dir ${CHRONO_CORE_REL_DLL} DIRECTORY)
        message(STATUS "  Chrono release DLL directory:     ${dir}")
        list(APPEND DLL_DIRS ${dir})
    endif()
    get_target_property(CHRONO_CORE_DBG_DLL Chrono::Chrono_core IMPORTED_LOCATION_DEBUG)
    if(CHRONO_CORE_DBG_DLL)
        get_filename_component(dir ${CHRONO_CORE_DBG_DLL} DIRECTORY)
        message(STATUS "  Chrono debug DLL directory:       ${dir}")
        list(APPEND DLL_DIRS ${dir})
    endif()

    # HDF5 DLL path
    if(TARGET hdf5::hdf5-static)
        message(STATUS "  No HDF5 DLL (hdf5-static)")
    elseif(TARGET hdf5::hdf5-shared)
        get_target_property(dll hdf5::hdf5-shared IMPORTED_LOCATION_RELEASE)
        get_filename_component(dir ${dll} DIRECTORY)
        message(STATUS "  HDF5 DLL directory (hdf5-shared): ${dir}")
        list(APPEND DLL_DIRS ${dir})
    elseif(TARGET HDF5::HDF5)
        get_target_property(dll HDF5::HFDF5 IMPORTED_LOCATION_RELEASE)
        get_filename_component(dir ${dll} DIRECTORY)
        message(STATUS "  HDF5 DLL directory (HDF5::HDF5):  ${dir}")
        list(APPEND DLL_DIRS ${dir})
    elseif(DEFINED HDF5_ROOT AND EXISTS "${HDF5_ROOT}/bin")
        message(STATUS "  HDF5 DLL directory (HDF5_ROOT):   ${HDF5_ROOT}/bin")
        list(APPEND DLL_DIRS ${HDF5_ROOT}/bin)
    else()
        message(STATUS "  Cannot infer HDF5 DLL directory")
    endif()

    # Irrlicht DLL path
    if(HYDROCHRONO_ENABLE_IRRLICHT)
        message(STATUS "  Irrlicht DLL directory:           ${Irrlicht_DLL_DIR}")
        list(APPEND DLL_DIRS ${Irrlicht_DLL_DIR})
    endif()

    # VSG DLL path
    if(HYDROCHRONO_ENABLE_VSG)
        message(STATUS "  VSG DLL directory:                ${VSG_DLL_DIR}")
        list(APPEND DLL_DIRS ${VSG_DLL_DIR})
    endif()

    set(TEST_ENVIRONMENT "PATH=${DLL_DIRS};$ENV{PATH}" PARENT_SCOPE)
endfunction()

# YAML-driven CLI app
if(HYDROCHRONO_ENABLE_YAML_RUNNER)
    add_subdirectory(app)
else()
    message(STATUS "YAML runner disabled - not creating run_hydrochrono target")
endif()

# Demos
if(HYDROCHRONO_ENABLE_DEMOS)
    add_subdirectory(demos)
endif()

# Tests
if(HYDROCHRONO_ENABLE_TESTS)
    # Use environment variable to override default data directory
    if(DEFINED ENV{HYDROCHRONO_DATA_DIR})
        set(HYDROCHRONO_DATA_DIR $ENV{HYDROCHRONO_DATA_DIR})
    else()
        set(HYDROCHRONO_DATA_DIR "${PROJECT_SOURCE_DIR}/demos")
    endif()

    enable_testing()
    
    # Configure test environment with DLL paths
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        configure_test_environment()
    endif()

    add_subdirectory(tests/regression)
    add_subdirectory(tests/unit)
endif()

# ===============================================================================
# ------ Installation -----------------------------------------------------------
# ===============================================================================

include(GNUInstallDirs)
option(HC_INSTALL_DEV_KIT "Install CMake dev package (headers/libs)" OFF)

if(HC_INSTALL_DEV_KIT)
    # Export HydroChrono as an importable target for other CMake projects
    install(TARGETS HydroChrono
        EXPORT HydroChronoTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(EXPORT HydroChronoTargets
        FILE HydroChronoTargets.cmake
        DESTINATION "${CMAKE_INSTALL_DATADIR}/HydroChrono/cmake"
        NAMESPACE HydroChrono::
    )

    # Install public headers
    install(DIRECTORY include/hydroc
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    
    # Install config.h header
    install(FILES ${PROJECT_BINARY_DIR}/hydroc/config.h-install
        RENAME "config.h"
        DESTINATION include/hydroc)

    # Install version.h header
    install(FILES ${PROJECT_BINARY_DIR}/hydroc/version.h
        DESTINATION include/hydroc)

    # Generate and install package config files for find_package() support
    include(CMakePackageConfigHelpers)
    configure_package_config_file(
        cmake/HydroChronoConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/HydroChronoConfig.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_DATADIR}/HydroChrono/cmake"
    )

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/HydroChronoConfig.cmake"
        DESTINATION "${CMAKE_INSTALL_DATADIR}/HydroChrono/cmake"
    )
endif()

if(HC_INSTALL_DEV_KIT)
    # Export targets for build tree usage (without installation)
    export(EXPORT HydroChronoTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/HydroChronoTargets.cmake"
        NAMESPACE HydroChrono::
    )
endif()

if(MSVC)
    # Windows: Enable debug symbols in Release builds with MSVC
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi /O2")
else()
    # Non-MSVC (Linux, macOS): Use -g for debug symbols
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O2")
endif()

# ------------------------------------------------------------------------------
# Runtime installer and ZIP packaging for releases
# ------------------------------------------------------------------------------

# Flat install tree for public distribution
option(HC_INSTALL_DEV_DEMOS "Install developer demo executables" OFF)

# Install main CLI and project shared libraries (if built)
if(TARGET run_hydrochrono)
    install(TARGETS run_hydrochrono  CONFIGURATIONS Release DESTINATION bin COMPONENT runtime)
endif()
if(TARGET HydroChrono)
    install(TARGETS HydroChrono      CONFIGURATIONS Release RUNTIME DESTINATION bin COMPONENT runtime)
endif()
if(TARGET HydroChronoGUI)
    install(TARGETS HydroChronoGUI   CONFIGURATIONS Release RUNTIME DESTINATION bin COMPONENT runtime)
endif()

# On Windows, install DLLs
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    message(STATUS "Set DLLs for installation")

    # Chrono DLLs — glob ALL DLLs from the Chrono bin directory to catch
    # transitive dependencies (e.g., Chrono_fsitdpf_vsg, yaml-cpp) that may
    # not be listed in CHRONO_TARGETS.
    get_target_property(_chrono_core_dll Chrono::Chrono_core IMPORTED_LOCATION_RELEASE)
    if(_chrono_core_dll)
        get_filename_component(_chrono_dll_dir "${_chrono_core_dll}" DIRECTORY)
        file(GLOB _chrono_dlls "${_chrono_dll_dir}/*.dll")
        foreach(_dll ${_chrono_dlls})
            message(STATUS "  Chrono DLL ${_dll}")
            install(FILES "${_dll}" DESTINATION bin COMPONENT runtime)
        endforeach()
    endif()

    # HDF5 DLLs
    foreach(tgt ${HDF5_TARGETS})
        get_target_property(tgt_DLL ${tgt} IMPORTED_LOCATION_RELEASE)
        if(EXISTS ${tgt_DLL})
            message(STATUS "  HDF5 DLL ${tgt_DLL}")
            install(FILES ${tgt_DLL} DESTINATION bin COMPONENT runtime)
        endif()
    endforeach()
    
    # Irrlicht DLL
    if(HYDROCHRONO_ENABLE_IRRLICHT AND EXISTS ${Irrlicht_DLL_DIR})
        set(_irr_dll "${Irrlicht_DLL_DIR}/Irrlicht.dll")
        if(EXISTS "${_irr_dll}")
            message(STATUS "  Irrlicht DLL ${_irr_dll}")
            install(FILES "${_irr_dll}" DESTINATION bin COMPONENT runtime)
        endif()
    endif()

    # VSG DLLs
    if(HYDROCHRONO_ENABLE_VSG AND EXISTS "${VSG_DLL_DIR}")
        file(GLOB _vsg_dlls "${VSG_DLL_DIR}/*.dll")
        if(_vsg_dlls)
            message(STATUS "  VSG DLLs ${_vsg_dlls}")
            install(FILES ${_vsg_dlls} DESTINATION bin COMPONENT runtime)
        endif()
    endif()
    
    # VSG DLLs (vsg, vsgImGui, vsgXchange, draco, etc.)
    if(HYDROCHRONO_ENABLE_VSG)
        # VSG_DLL_DIR is set by Chrono's find_package when VSG is found
        # Fall back to deriving from vsg_DIR if VSG_DLL_DIR is not set
        if(NOT DEFINED VSG_DLL_DIR OR VSG_DLL_DIR STREQUAL "")
            if(DEFINED vsg_DIR)
                # vsg_DIR points to lib/cmake/vsg, go up 3 levels to root then into bin
                get_filename_component(_vsg_root "${vsg_DIR}" DIRECTORY)  # lib/cmake
                get_filename_component(_vsg_root "${_vsg_root}" DIRECTORY)  # lib
                get_filename_component(_vsg_root "${_vsg_root}" DIRECTORY)  # root
                set(VSG_DLL_DIR "${_vsg_root}/bin")
            endif()
        endif()
        
        if(DEFINED VSG_DLL_DIR AND EXISTS "${VSG_DLL_DIR}")
            file(GLOB _vsg_dlls "${VSG_DLL_DIR}/*.dll")
            if(_vsg_dlls)
                message(STATUS "Installing VSG DLLs from: ${VSG_DLL_DIR}")
                install(FILES ${_vsg_dlls} DESTINATION bin COMPONENT runtime)
            endif()
        else()
            message(WARNING "HYDROCHRONO_ENABLE_VSG is ON but VSG_DLL_DIR not found. VSG DLLs will not be packaged.")
        endif()
    endif()
    
    # Python DLLs, if Chrono::Parsers depends on it.
    # Derive the DLL path from the already-found interpreter rather than calling
    # find_package(Python3 ... Development), which can stall on Windows/conda.
    if(CHRONO_PARSERS_PYTHON AND Python3_EXECUTABLE)
        get_filename_component(_py_bin_dir "${Python3_EXECUTABLE}" DIRECTORY)
        file(GLOB _py_dlls "${_py_bin_dir}/python3*.dll")
        foreach(_py_dll ${_py_dlls})
            message(STATUS "  Python DLL ${_py_dll}")
            install(FILES "${_py_dll}" DESTINATION bin COMPONENT runtime)
        endforeach()
    endif()

endif()

# Public, user-facing regression test suite only
install(DIRECTORY ${PROJECT_SOURCE_DIR}/tests/regression/run_hydrochrono
        DESTINATION tests COMPONENT python-tests
        PATTERN "__pycache__" EXCLUDE
        PATTERN "*.pyc" EXCLUDE)

# Install test input data (YAML configs, hydro data, geometry) into tests/run_hydrochrono
# The Python scripts expect test case folders (iea_sphere/, oswec/, rm3/, f3of/) to be siblings
install(DIRECTORY ${PROJECT_SOURCE_DIR}/data/demos/run_hydrochrono/
        DESTINATION tests/run_hydrochrono COMPONENT python-tests
        PATTERN "__pycache__" EXCLUDE
        PATTERN "*.pyc" EXCLUDE)

# Simple runner script to execute tests from the installed tree
install(PROGRAMS ${PROJECT_SOURCE_DIR}/scripts/RUN-TESTS.ps1
        DESTINATION tests COMPONENT python-tests)

# MSVC runtime DLLs and ZIP packaging via CPack
# InstallRequiredSystemLibraries handles msvcp/vcruntime/ucrt, but not vcomp (OpenMP).
# Append the OpenMP runtime so it's included in the package.
if(MSVC AND OpenMP_CXX_FOUND)
    get_filename_component(_msvc_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
    find_file(_vcomp_dll "vcomp140.dll"
        PATHS "${_msvc_dir}" "${_msvc_dir}/../../../redist/x64" "C:/Windows/System32"
        NO_DEFAULT_PATH)
    if(_vcomp_dll)
        list(APPEND CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS "${_vcomp_dll}")
        message(STATUS "  OpenMP DLL ${_vcomp_dll}")
    else()
        message(WARNING "vcomp140.dll not found — OpenMP runtime will be missing from package")
    endif()
endif()
include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "ZIP")
set(CPACK_PACKAGE_NAME "HydroChrono")
set(CPACK_COMPONENTS_ALL runtime python-tests dev-demos)
include(CPack)
