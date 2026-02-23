/*********************************************************************
 * @file  setup_from_yaml.cpp
 * @brief Creates HydroSystem from parsed YAML configuration.
 *
 * Implements SetupHydroFromYAML(): matches YAML body configs to Chrono
 * bodies, creates wave models, and returns a configured HydroSystem.
 *********************************************************************/

#include "setup_from_yaml.h"
#include "config_loader.h"
#include <hydroc/hydro_system.h> // For HydroSystem
#ifdef HYDROCHRONO_HAVE_MOORDYN
#include <hydroc/config/moordyn_config.h>
#endif
#include <hydroc/waves/wave_base.h>
#include <hydroc/waves/regular_wave.h>
#include <hydroc/waves/irregular_wave.h>
#include <hydroc/logging.h>         // For Logger
#include "../radiation/radiation_rirf_processing.h" // For TaperedDirectOptions (canonical type)
#include "../force_components/radiation_component.h" // For RadiationConvolutionMode (canonical type)
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <cmath> // For M_PI
#include <unordered_map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace chrono;

namespace {

// ------------------------------------------------------------
// SECTION: Wave model factory
// ------------------------------------------------------------
// Creates appropriate WaveBase implementation from YAML settings.
// TODO: Extract to dedicated wave factory module in future refactor.

/**
 * @brief Create a wave object from wave settings.
 */
std::shared_ptr<WaveBase> CreateWaveFromSettings(const WaveSettings& wave_settings,
                                                 double timestep,
                                                 double sim_duration,
                                                 double ramp_duration) {
    // Normalize type to lowercase to be tolerant of input casing
    std::string type = wave_settings.type;
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    if (type == "regular") {
        auto regular_wave = std::make_shared<RegularWave>();
        regular_wave->SetAmplitude(wave_settings.height / 2.0);
        regular_wave->SetPeriod(wave_settings.period);
        regular_wave->SetPhase(wave_settings.phase);
        
        hydroc::debug::LogDebug(std::string("Attached wave model: RegularWave, H=") + std::to_string(wave_settings.height) + 
                            "m, T=" + std::to_string(wave_settings.period) + "s");
        
        return regular_wave;
        
    } else if (type == "irregular") {
        IrregularWaveParams params;
        params.ramp_duration = ramp_duration;
        params.wave_height = wave_settings.height;
        params.wave_period = wave_settings.period;
        params.nfrequencies = 1000;
        params.seed = (wave_settings.seed > 0 ? wave_settings.seed : 1);
        
        auto irregular_wave = std::make_shared<IrregularWaves>(params);
        
        hydroc::debug::LogDebug(std::string("Attached wave model: IrregularWaves, H=") + std::to_string(wave_settings.height) + 
                            "m, T=" + std::to_string(wave_settings.period) + "s");
        
        return irregular_wave;
        
    } else if (type == "no_wave" || type == "still_ci" || type == "still") {
        auto no_wave = std::make_shared<NoWave>();
        
        hydroc::debug::LogDebug("Attached wave model: NoWave (still water)");
        
        return no_wave;
        
    } else {
        throw std::runtime_error("Unsupported wave type: " + wave_settings.type);
    }
}

/**
 * @brief Match hydrodynamic bodies with Chrono bodies by name.
 */
std::vector<std::shared_ptr<ChBody>> MatchBodiesByName(
    const std::vector<HydroBody>& hydro_bodies,
    const std::vector<std::shared_ptr<ChBody>>& chrono_bodies,
    std::string& h5_file_path) {
    
    std::vector<std::shared_ptr<ChBody>> matched_bodies;
    
    // For now, we'll use the first H5 file found (assuming all bodies use the same file)
    // In the future, this could be enhanced to support different H5 files per body
    if (!hydro_bodies.empty()) {
        h5_file_path = hydro_bodies[0].h5_file;
    }
    
    // Match bodies by name
    for (const auto& hydro_body : hydro_bodies) {
        bool found = false;
        
        for (const auto& chrono_body : chrono_bodies) {
            if (chrono_body->GetName() == hydro_body.name) {
                matched_bodies.push_back(chrono_body);
                found = true;
                
                // Log the matched body details
                hydroc::debug::LogDebug(std::string("Body: ") + hydro_body.name + 
                          " -> h5: " + h5_file_path + 
                          ", excitation: " + (hydro_body.include_excitation ? "true" : "false") + 
                          ", radiation: " + (hydro_body.include_radiation ? "true" : "false"));
                
                break;
            }
        }
        
        if (!found) {
            hydroc::cli::LogWarning("Hydrodynamic body '" + hydro_body.name + "' not found in Chrono system");
        }
    }
    
    return matched_bodies;
}

// ------------------------------------------------------------
// SECTION: Body matching
// ------------------------------------------------------------
// Matches YAML body configurations to Chrono bodies by name.
// TODO: Extract to separate module, support per-body H5 files.

} // anonymous namespace

// ------------------------------------------------------------
// SECTION: Main setup function
// ------------------------------------------------------------
// Orchestrates body matching, wave creation, and HydroSystem initialization.

std::unique_ptr<HydroSystem> SetupHydroFromYAML(
    const YAMLHydroData& hydro_data,
    const std::vector<std::shared_ptr<ChBody>>& bodies,
    double timestep,
    double sim_duration,
    double ramp_duration) {
    
    // Match hydrodynamic bodies with Chrono bodies (multibody: matches all configured bodies)
    std::string h5_file_path;
    auto matched_bodies = MatchBodiesByName(hydro_data.bodies, bodies, h5_file_path);
    
    if (matched_bodies.empty()) {
        throw std::runtime_error("No hydrodynamic bodies found in Chrono system");
    }
    
    // Create wave object from settings (system-wide, not per-body)
    auto wave = CreateWaveFromSettings(hydro_data.waves, timestep, sim_duration, ramp_duration);
    
    // Create and initialize HydroSystem (multibody: all matched bodies passed in)
    auto hydro_system = std::make_unique<HydroSystem>(matched_bodies, h5_file_path, wave);
    
    hydroc::debug::LogDebug(std::string("Initialized HydroSystem with ") + std::to_string(matched_bodies.size()) + " bodies");

    // ─────────────────────────────────────────────────────────────────────────
    // Radiation method selection
    // ─────────────────────────────────────────────────────────────────────────
    std::string method = hydro_data.radiation_method;
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);
    
    if (method == "state_space") {
        hydro_system->SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kStateSpace);
        hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Radiation Method", "StateSpace"));
        
        // Set state-space options
        hydrochrono::hydro::StateSpaceOptions ss_opts;
        ss_opts.max_order = hydro_data.ss_max_order;
        ss_opts.r2_threshold = hydro_data.ss_r2_threshold;
        ss_opts.max_hankel_size = hydro_data.ss_max_hankel_size;
        ss_opts.r2_num_samples = hydro_data.ss_r2_num_samples;
        hydro_system->SetStateSpaceOptions(ss_opts);
        
        // Enable kernel fit diagnostics if requested
        if (hydro_data.output_kernel_fit) {
            hydro_system->SetOutputKernelFit(true);
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Kernel Fit Diagnostics", "Enabled"));
        }
        
        if (hydroc::debug::IsDebugEnabled()) {
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "SS Max Order", std::to_string(ss_opts.max_order)));
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "SS R² Threshold", std::to_string(ss_opts.r2_threshold)));
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "SS Max Hankel Size", std::to_string(ss_opts.max_hankel_size)));
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "SS R² Samples", std::to_string(ss_opts.r2_num_samples)));
        }
    } else {
        // Default: rirf_convolution (or any unrecognized value)
        hydro_system->SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kRirfConvolution);
        hydroc::debug::LogDebug("Radiation method: RIRF Convolution (default)");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Convolution settings (applies when method == rirf_convolution)
    // ─────────────────────────────────────────────────────────────────────────
    std::string mode = hydro_data.radiation_convolution_mode;
    hydroc::debug::LogDebug("Parsed convolution mode: '" + mode + "'");
    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
    hydroc::debug::LogDebug("Lowercase mode: '" + mode + "'");
    if (mode == "tapereddirect") {
        hydro_system->SetRadiationConvolutionMode(hydrochrono::hydro::RadiationConvolutionMode::TaperedDirect);
        hydroc::debug::LogDebug("Radiation convolution mode: TaperedDirect");
        hydrochrono::hydro::TaperedDirectOptions opts;
        opts.smoothing = !hydro_data.td_smoothing.empty() ? hydro_data.td_smoothing : opts.smoothing;
        opts.window_length = std::max(3, hydro_data.td_window_length != 0 ? hydro_data.td_window_length : opts.window_length);
        if (opts.window_length % 2 == 0) opts.window_length += 1; // enforce odd
        
        // RIRF truncation
        opts.rirf_end_time = hydro_data.td_rirf_end_time;
        
        // Simple taper control
        opts.taper_start_percent = hydro_data.td_taper_start_percent;
        opts.taper_end_percent = hydro_data.td_taper_end_percent;
        opts.taper_final_amplitude = hydro_data.td_taper_final_amplitude;
        opts.export_plot_csv = hydro_data.td_export_plot_csv;
        hydro_system->SetTaperedDirectOptions(opts);

        // CLI inline bullets near main summary
        hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Convolution Mode", "TaperedDirect"));
        if (hydroc::debug::IsDebugEnabled()) {
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Conv Smoothing", opts.smoothing));
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Conv Window Length", std::to_string(opts.window_length)));
            if (opts.rirf_end_time > 0.0) {
                hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Conv RIRF End Time", std::to_string(opts.rirf_end_time) + "s"));
            }
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Conv Taper Start %", std::to_string(opts.taper_start_percent)));
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Conv Taper End %", std::to_string(opts.taper_end_percent)));
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Conv Taper Final Amp", std::to_string(opts.taper_final_amplitude)));
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Conv Export CSV", (opts.export_plot_csv ? "true" : "false")));
        }
    } else {
        hydro_system->SetRadiationConvolutionMode(hydrochrono::hydro::RadiationConvolutionMode::Baseline);
        hydroc::debug::LogDebug("Radiation convolution mode: Baseline");
        hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Convolution Mode", "Baseline"));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // MoorDyn mooring coupling (optional)
    // ─────────────────────────────────────────────────────────────────────────
#ifdef HYDROCHRONO_HAVE_MOORDYN
    if (hydro_data.moordyn_enabled && !hydro_data.moordyn_input_file.empty()) {
        MoorDynConfig md_cfg;
        md_cfg.enabled = true;
        md_cfg.input_file = hydro_data.moordyn_input_file;

        // Resolve body names to 0-based indices in matched_bodies
        for (const auto& bname : hydro_data.moordyn_body_names) {
            bool found = false;
            for (size_t i = 0; i < matched_bodies.size(); ++i) {
                if (matched_bodies[i]->GetName() == bname) {
                    md_cfg.coupled_body_indices.push_back(static_cast<int>(i));
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error(
                    "MoorDyn config references body '" + bname +
                    "' which was not found among hydrodynamic bodies");
            }
        }

        hydro_system->SetMoorDynConfig(md_cfg);
        hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine(
            "+", "MoorDyn", hydro_data.moordyn_input_file +
            " (" + std::to_string(md_cfg.coupled_body_indices.size()) + " bodies)"));
    }
#else
    if (hydro_data.moordyn_enabled) {
        hydroc::cli::LogWarning(
            "MoorDyn is enabled in YAML but HydroChrono was built without "
            "HYDROCHRONO_ENABLE_MOORDYN. Mooring forces will be ignored.");
    }
#endif

    return hydro_system;
}

// ------------------------------------------------------------
// SECTION: Convenience loader + setup helper
// ------------------------------------------------------------
// Keeps file parsing and setup logic together for callers that
// only know the YAML path.

std::unique_ptr<HydroSystem> SetupHydroFromYAMLFile(
    const std::string& hydro_yaml_path,
    const std::vector<std::shared_ptr<ChBody>>& bodies,
    double timestep,
    double sim_duration,
    double ramp_duration) {
    const YAMLHydroData hydro_data = hydrochrono::hydro::LoadHydroConfigFromYaml(hydro_yaml_path);
    return SetupHydroFromYAML(hydro_data, bodies, timestep, sim_duration, ramp_duration);
}
