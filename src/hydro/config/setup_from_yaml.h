/*********************************************************************
 * @file  setup_from_yaml.h
 * @brief Factory function to create HydroForces from YAML configuration.
 *
 * Provides SetupHydroFromYAML() which matches YAML body configs to Chrono
 * bodies, creates wave models, and returns a configured HydroForces.
 *********************************************************************/

#ifndef HYDRO_CONFIG_SETUP_FROM_YAML_H
#define HYDRO_CONFIG_SETUP_FROM_YAML_H

#include <hydroc/config/hydro_config.h>
#include <chrono/physics/ChBody.h>
#include <memory>
#include <vector>
#include <string>

// Forward declarations
class HydroForces;

/**
 * @brief Setup hydrodynamic forces from parsed YAML data.
 *
 * This function connects parsed YAML data (from hydro.yaml) to actual hydrodynamic forces.
 * It builds the appropriate WaveBase subclass from hydro_data.waves, matches body names
 * with their corresponding HDF5 files, and initializes HydroForces with the matched bodies.
 *
 * @param hydro_data Parsed hydrodynamic configuration from hydro.yaml
 * @param bodies Vector of Chrono bodies that may have hydrodynamic forces
 * @param timestep Simulation timestep (used for irregular wave setup)
 * @param sim_duration Simulation duration (used for irregular wave setup)
 * @param ramp_duration Wave ramp duration (used for irregular wave setup)
 * @return Unique pointer to initialized HydroForces object
 * @throws std::runtime_error on configuration errors
 */
std::unique_ptr<HydroForces> SetupHydroFromYAML(
    const YAMLHydroData& hydro_data,
    const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
    double timestep,
    double sim_duration,
    double ramp_duration
);

/**
 * @brief Convenience overload that loads the YAML file and then sets up HydroForces.
 *
 * Behaviour matches calling LoadHydroConfigFromYaml followed by SetupHydroFromYAML.
 *
 * @param hydro_yaml_path Path to hydro.yaml.
 * @param bodies Chrono bodies that may receive hydrodynamic forces.
 * @param timestep Simulation time step (used for irregular waves).
 * @param sim_duration Simulation duration (irregular waves).
 * @param ramp_duration Wave ramp duration.
 * @return Initialized HydroForces object.
 */
std::unique_ptr<HydroForces> SetupHydroFromYAMLFile(
    const std::string& hydro_yaml_path,
    const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
    double timestep,
    double sim_duration,
    double ramp_duration
);

#endif  // HYDRO_CONFIG_SETUP_FROM_YAML_H

