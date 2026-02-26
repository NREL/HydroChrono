/*********************************************************************
 * @file  setup_from_yaml.cpp
 * @brief Creates HydroSystem from parsed YAML configuration.
 *
 * Implements SetupHydroFromYAML(): matches YAML body configs to Chrono
 * bodies, creates wave models, and returns a configured HydroSystem.
 *********************************************************************/

#include "setup_from_yaml.h"
#include "config_loader.h"
#include <hydroc/hydro_system.h>
#ifdef HYDROCHRONO_HAVE_MOORDYN
#include <hydroc/config/moordyn_config.h>
#endif
#include <hydroc/waves/wave_base.h>
#include <hydroc/waves/regular_wave.h>
#include <hydroc/waves/irregular_wave.h>
#include <hydroc/logging.h>
#include <hydroc/radiation/radiation_types.h>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

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
        params.ramp_duration = (wave_settings.ramp_duration > 0.0) ? wave_settings.ramp_duration : ramp_duration;
        params.wave_height   = wave_settings.height;
        params.wave_period   = wave_settings.period;
        params.nfrequencies  = (wave_settings.nfrequencies > 0)
                               ? wave_settings.nfrequencies
                               : IrregularWaveParams::kDefaultNFrequencies;
        params.seed          = (wave_settings.seed > 0 ? wave_settings.seed : 1);
        params.eta_file_path = wave_settings.eta_file;
        if (wave_settings.frequency_min > 0.0) params.frequency_min = wave_settings.frequency_min;
        if (wave_settings.frequency_max > 0.0) params.frequency_max = wave_settings.frequency_max;

        std::string ramp_lower = wave_settings.ramp_type;
        std::transform(ramp_lower.begin(), ramp_lower.end(), ramp_lower.begin(), ::tolower);
        if (ramp_lower == "cosine")
            params.ramp_type = ExcitationRampType::kCosine;

        auto irregular_wave = std::make_shared<IrregularWaves>(params);

        if (!wave_settings.eta_file.empty()) {
            hydroc::debug::LogDebug(std::string("Attached wave model: IrregularWaves (elevation import), eta_file=") +
                                wave_settings.eta_file);
        } else {
            hydroc::debug::LogDebug(std::string("Attached wave model: IrregularWaves, H=") + std::to_string(wave_settings.height) +
                                "m, T=" + std::to_string(wave_settings.period) + "s");
        }

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
    // Convolution truncation (excitation + radiation)
    // ─────────────────────────────────────────────────────────────────────────
    if (hydro_data.excitation_truncation_time > 0.0) {
        hydro_system->SetExcitationTruncationTime(hydro_data.excitation_truncation_time);
        hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine(
            "•", "Excitation Truncation", std::to_string(hydro_data.excitation_truncation_time) + "s"));
    }
    if (hydro_data.radiation_truncation_time > 0.0) {
        hydro_system->SetRadiationTruncationTime(hydro_data.radiation_truncation_time);
        hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine(
            "•", "Radiation Truncation", std::to_string(hydro_data.radiation_truncation_time) + "s"));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Radiation method selection
    // ─────────────────────────────────────────────────────────────────────────
    std::string method = hydro_data.radiation_method;
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);
    
    if (method == "state_space") {
        hydro_system->SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kStateSpace);
        hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Radiation Method", "StateSpace"));
        
        hydrochrono::hydro::StateSpaceOptions ss_opts;
        ss_opts.max_order = hydro_data.ss_max_order;
        ss_opts.r2_threshold = hydro_data.ss_r2_threshold;
        ss_opts.max_hankel_size = hydro_data.ss_max_hankel_size;
        ss_opts.r2_num_samples = hydro_data.ss_r2_num_samples;
        hydro_system->SetStateSpaceOptions(ss_opts);
        
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
        hydro_system->SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kRirfConvolution);
        hydroc::debug::LogDebug("Radiation method: RIRF Convolution (default)");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Radiation kernel processing (smoothing + tapering, composable)
    // ─────────────────────────────────────────────────────────────────────────
    auto kp = hydro_data.radiation_kernel_processing;
    kp.smoothing_window = std::max(3, kp.smoothing_window);
    if (kp.smoothing_window % 2 == 0) kp.smoothing_window += 1;  // enforce odd

    if (kp.RequiresProcessing()) {
        hydro_system->SetRadiationKernelProcessing(kp);
        if (kp.smoothing_type != "none") {
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "RIRF Smoothing", kp.smoothing_type));
        }
        if (kp.taper_enabled) {
            hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "RIRF Taper", "enabled"));
            if (hydroc::debug::IsDebugEnabled()) {
                hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Taper Start %", std::to_string(kp.taper_start_percent)));
                hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Taper End %", std::to_string(kp.taper_end_percent)));
                hydroc::cli::LogInfo(hydroc::cli::CreateAlignedLine("•", "Taper Final Amp", std::to_string(kp.taper_final_amplitude)));
            }
        }
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
