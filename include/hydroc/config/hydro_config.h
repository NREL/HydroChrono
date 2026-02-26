/*********************************************************************
 * @file  hydro_config.h
 * @brief Public configuration types for parsed YAML hydro settings.
 *
 * MAIN TYPES:
 *   - HydroBody: Per-body hydrodynamic settings (H5 file, flags)
 *   - WaveSettings: Wave type, height, period, spectrum
 *   - YAMLHydroData: Top-level container for hydro.yaml data
 *********************************************************************/

#ifndef HYDROC_CONFIG_HYDRO_CONFIG_H
#define HYDROC_CONFIG_HYDRO_CONFIG_H

#include <hydroc/radiation/radiation_types.h>
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
    // Elevation import (if non-empty, overrides spectral generation)
    std::string eta_file;
    // Excitation ramp shape: "linear" (default) or "cosine" (WEC-Sim convention)
    std::string ramp_type    = "linear";
    double ramp_duration     = 0.0;   // Ramp duration [s]; 0 = no ramp (overrides function param if > 0)
    double frequency_min  = 0.0;   // 0 = use IrregularWaveParams::kDefaultFreqMin
    double frequency_max  = 0.0;   // 0 = use IrregularWaveParams::kDefaultFreqMax
    int    nfrequencies   = 0;     // 0 = use IrregularWaveParams::kDefaultNFrequencies
    // TODO: Add spectrum parameters (peak enhancement factor, etc.)
};

/**
 * @brief Top-level container for hydrodynamic configuration data from YAML.
 */
struct YAMLHydroData {
    std::vector<HydroBody> bodies;
    WaveSettings waves;
    
    // ─────────────────────────────────────────────────────────────────────────
    // Convolution truncation (system-wide, wave-type-independent)
    // ─────────────────────────────────────────────────────────────────────────
    double excitation_truncation_time = 0.0;  ///< Excitation IRF truncation [s]; 0 = full
    double radiation_truncation_time = 0.0;   ///< Radiation RIRF truncation [s]; 0 = full
    
    // ─────────────────────────────────────────────────────────────────────────
    // Radiation method selection (system-wide)
    // ─────────────────────────────────────────────────────────────────────────
    std::string radiation_method = "rirf_convolution";  ///< "rirf_convolution" or "state_space"
    
    // ─────────────────────────────────────────────────────────────────────────
    // Radiation kernel processing (smoothing + tapering, composable)
    // ─────────────────────────────────────────────────────────────────────────
    hydrochrono::hydro::RadiationKernelProcessing radiation_kernel_processing;
    
    // ─────────────────────────────────────────────────────────────────────────
    // State-space options (only used if radiation_method == "state_space")
    // ─────────────────────────────────────────────────────────────────────────
    int ss_max_order = 10;
    double ss_r2_threshold = 0.95;
    int ss_max_hankel_size = 200;
    int ss_r2_num_samples = 50;
    
    // ─────────────────────────────────────────────────────────────────────────
    // Diagnostics
    // ─────────────────────────────────────────────────────────────────────────
    bool output_kernel_fit = false;

    // ─────────────────────────────────────────────────────────────────────────
    // MoorDyn mooring coupling (optional)
    // ─────────────────────────────────────────────────────────────────────────
    bool moordyn_enabled = false;
    std::string moordyn_input_file;
    std::vector<std::string> moordyn_body_names;  // e.g. ["body1"]
};

#endif  // HYDROC_CONFIG_HYDRO_CONFIG_H
