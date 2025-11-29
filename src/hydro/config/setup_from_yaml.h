/*********************************************************************
 * @file  setup_from_yaml.h
 *
 * @brief Setup hydrodynamic forces from parsed YAML data.
 *
 * OVERVIEW:
 * This file provides the main entry point for configuring hydrodynamic
 * forces from YAML configuration files. It bridges parsed configuration
 * data to the TestHydro force computation system.
 *
 * MAIN RESPONSIBILITIES:
 * - Match Chrono bodies with YAML body configurations by name
 * - Create appropriate WaveBase implementations from wave settings
 * - Initialize TestHydro with matched bodies and wave models
 * - Configure radiation convolution options (Baseline/TaperedDirect)
 *
 * INTERACTIONS:
 * - Takes parsed YAMLHydroData from yaml_parser
 * - Takes Chrono bodies from the simulation system
 * - Returns configured TestHydro instance ready for simulation
 * - Uses wave_types.h for wave model creation
 *
 * KEY ASSUMPTIONS:
 * - All bodies share the same H5 file (uses first body's h5_file)
 * - Body names in YAML match Chrono body names exactly
 * - Bodies are matched in order found in YAML (not sorted)
 * - Wave settings apply system-wide (not per-body)
 *
 * KNOWN LIMITATIONS:
 * - Single H5 file assumption (no per-body H5 files)
 * - No validation of body count vs H5 file contents
 * - Wave configuration is system-wide only
 * - Body matching is simple string equality (no fuzzy matching)
 *********************************************************************/

#ifndef HYDRO_CONFIG_SETUP_FROM_YAML_H
#define HYDRO_CONFIG_SETUP_FROM_YAML_H

#include <hydroc/config/hydro_config.h>
#include <chrono/physics/ChBody.h>
#include <memory>
#include <vector>
#include <string>

// Forward declarations
class TestHydro;

/**
 * @brief Setup hydrodynamic forces from parsed YAML data.
 *
 * This function connects parsed YAML data (from hydro.yaml) to actual hydrodynamic forces.
 * It builds the appropriate WaveBase subclass from hydro_data.waves, matches body names
 * with their corresponding HDF5 files, and initializes TestHydro with the matched bodies.
 *
 * @param hydro_data Parsed hydrodynamic configuration from hydro.yaml
 * @param bodies Vector of Chrono bodies that may have hydrodynamic forces
 * @param timestep Simulation timestep (used for irregular wave setup)
 * @param sim_duration Simulation duration (used for irregular wave setup)
 * @param ramp_duration Wave ramp duration (used for irregular wave setup)
 * @return Unique pointer to initialized TestHydro object
 * @throws std::runtime_error on configuration errors
 */
std::unique_ptr<TestHydro> SetupHydroFromYAML(
    const YAMLHydroData& hydro_data,
    const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
    double timestep,
    double sim_duration,
    double ramp_duration
);

/**
 * @brief Convenience overload that loads the YAML file and then sets up TestHydro.
 *
 * Behaviour matches calling LoadHydroConfigFromYaml followed by SetupHydroFromYAML.
 *
 * @param hydro_yaml_path Path to hydro.yaml.
 * @param bodies Chrono bodies that may receive hydrodynamic forces.
 * @param timestep Simulation time step (used for irregular waves).
 * @param sim_duration Simulation duration (irregular waves).
 * @param ramp_duration Wave ramp duration.
 * @return Initialized TestHydro object.
 */
std::unique_ptr<TestHydro> SetupHydroFromYAMLFile(
    const std::string& hydro_yaml_path,
    const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
    double timestep,
    double sim_duration,
    double ramp_duration
);

#endif  // HYDRO_CONFIG_SETUP_FROM_YAML_H

