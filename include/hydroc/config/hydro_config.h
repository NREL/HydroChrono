/*********************************************************************
 * @file  hydro_config.h
 *
 * @brief Public configuration types for HydroChrono YAML configuration.
 *
 * This header defines the data structures used to represent parsed
 * hydro.yaml configuration files. Users and library consumers should
 * include this header to work with hydrodynamic configuration data.
 *
 * MAIN TYPES:
 * - HydroBody: Configuration for a single hydrodynamic body
 * - WaveSettings: Configuration for wave conditions
 * - YAMLHydroData: Top-level container for all configuration
 *
 * USAGE:
 * Include <hydroc/config/hydro_config.h> to access these types.
 * Use hydrochrono::hydro::LoadHydroConfigFromYaml() to parse YAML files.
 *********************************************************************/

#ifndef HYDROC_CONFIG_HYDRO_CONFIG_H
#define HYDROC_CONFIG_HYDRO_CONFIG_H

#include <string>
#include <vector>
#include <array>

/**
 * @brief Configuration for a hydrodynamic body.
 */
struct HydroBody {
    std::string name = "";
    std::string h5_file = "";
    bool include_excitation = true;
    bool include_radiation = true;
    std::string radiation_calculation = "convolution";  // "convolution" or "state_space"
    // Optional convolution mode for radiation kernel preprocessing: "Baseline" (default) or "TaperedDirect"
    std::string radiation_convolution_mode = "Baseline";
    // Optional TaperedDirect tuning
    std::string td_smoothing = "sg";            // "sg" or "moving_average"
    int td_window_length = 5;                    // odd >= 3
    double td_rms_threshold_factor = 0.02;       // e.g. 0.02
    double td_taper_fraction_remaining = 0.25;   // e.g. 0.25
    bool td_export_plot_csv = false;             // export CSV for before/after
    // TODO: Add nonlinear buoyancy fields
    // TODO: Add drag coefficient fields
};

/**
 * @brief Configuration for wave settings.
 */
struct WaveSettings {
    std::string type = "regular";  // "regular", "irregular", "no_wave"
    double height = 0.0;
    double period = 0.0;
    double direction = 0.0;  // degrees, 0 = positive x
    double phase = 0.0;
    std::string spectrum = "pierson_moskowitz";  // "pierson_moskowitz", "jonswap", etc.
    int seed = -1; // optional irregular seed; -1 means unset
    // Sweep support (expanded values) for period; if empty, use 'period'
    std::vector<double> period_values;
    // TODO: Add spectrum parameters (peak enhancement factor, etc.)
};

/**
 * @brief Top-level container for hydrodynamic configuration data from YAML.
 */
struct YAMLHydroData {
    std::vector<HydroBody> bodies;
    WaveSettings waves;
    // Optional system-wide convolution settings
    std::string radiation_convolution_mode = "Baseline"; // Baseline | TaperedDirect
    std::string td_smoothing = "sg";
    int td_window_length = 5;
    
    // RIRF truncation
    double td_rirf_end_time = -1.0;
    
    // Simple taper control - sensible defaults for improved stability
    double td_taper_start_percent = 0.8;      // start taper at 80% (taper last 20%)
    double td_taper_end_percent = 1.0;        // end taper at 100% of total time series
    double td_taper_final_amplitude = 0.0;    // final amplitude as fraction of original (0.0 = zero, 1.0 = no change)
    bool td_export_plot_csv = false;          // dump before/after CSV summaries (false by default)
};

#endif  // HYDROC_CONFIG_HYDRO_CONFIG_H

