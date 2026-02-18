/*********************************************************************
 * @file  irregular_wave.cpp
 * @brief Irregular wave model implementation.
 *********************************************************************/

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <hydroc/waves/irregular_wave.h>
#include "wave_utilities.h"

#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <unsupported/Eigen/Splines>

#include <hydroc/logging.h>

IrregularWaves::IrregularWaves(const IrregularWaveParams& params) : params_(params) {
    wave_stretching_ = params.wave_stretching_;
}

void IrregularWaves::InitializeIRFVectors() {
    ex_irf_sampled_.resize(num_bodies_);
    ex_irf_time_sampled_.resize(num_bodies_);
    ex_irf_width_sampled_.resize(num_bodies_);

    for (unsigned int b = 0; b < num_bodies_; b++) {
        ex_irf_sampled_[b]      = GetExcitationIRF(b);
        ex_irf_time_sampled_[b] = wave_info_[b].excitation_irf_time;
        CalculateWidthIRF();
    }

    if (params_.simulation_dt_ > 0.0) {
        ResampleIRF(params_.simulation_dt_);
    }

    if (!params_.eta_file_path_.empty()) {
        ReadEtaFromFile();
        
        // Check if eta file was read successfully
        if (time_data_.empty() || free_surface_elevation_sampled_.empty()) {
            std::cerr << "Error: Failed to read eta data from file: " << params_.eta_file_path_ << std::endl;
            spectrumCreated_ = false;
            return;
        }
        
        // Populate free_surface_time_sampled_ from time_data_, extending backwards
        // to cover negative times needed for IRF convolution (same as CreateFreeSurfaceElevation does)
        double t_irf_max = 0.0;
        for (size_t b = 0; b < ex_irf_time_sampled_.size(); b++) {
            if (ex_irf_time_sampled_[b].size() > 0) {
                double t_end = ex_irf_time_sampled_[b][ex_irf_time_sampled_[b].size() - 1];
                if (t_end > t_irf_max) t_irf_max = t_end;
            }
        }
        
        double eta_tmin = time_data_.front();
        double required_tmin = -t_irf_max;
        
        if (required_tmin < eta_tmin && time_data_.size() > 1) {
            // Extend backwards with zero elevation padding
            double dt = time_data_[1] - time_data_[0];
            int num_pad = static_cast<int>(std::ceil((eta_tmin - required_tmin) / dt));
            
            std::vector<double> extended_time(num_pad + time_data_.size());
            std::vector<double> extended_eta(num_pad + free_surface_elevation_sampled_.size());
            
            for (int i = 0; i < num_pad; ++i) {
                extended_time[i] = eta_tmin - (num_pad - i) * dt;
                extended_eta[i] = 0.0;
            }
            for (size_t i = 0; i < time_data_.size(); ++i) {
                extended_time[num_pad + i] = time_data_[i];
                extended_eta[num_pad + i] = free_surface_elevation_sampled_[i];
            }
            
            free_surface_time_sampled_ = std::move(extended_time);
            free_surface_elevation_sampled_ = std::move(extended_eta);
        } else {
            free_surface_time_sampled_ = time_data_;
        }
        
        spectrumCreated_ = false;
    } else if (params_.wave_height_ != 0.0 && params_.wave_period_ != 0.0) {
        CreateSpectrum();
        CreateFreeSurfaceElevation();
        spectrumCreated_ = true;
    }
}

std::vector<double> IrregularWaves::GetSpectrum() {
    if (!spectrumCreated_) {
        throw std::runtime_error("Spectrum has not been created. Initialize with wave height and period to create spectrum.");
    }
    return spectrum_;
}

std::vector<double> IrregularWaves::GetFreeSurfaceElevation() {
    return free_surface_elevation_sampled_;
}

std::vector<double> IrregularWaves::GetFreeSurfaceTime() const {
    return free_surface_time_sampled_;
}

std::vector<double> IrregularWaves::GetFrequenciesHz() const {
    std::vector<double> out(spectrum_frequencies_.size());
    for (Eigen::Index i = 0; i < spectrum_frequencies_.size(); ++i) {
        out[i] = spectrum_frequencies_[i];
    }
    return out;
}

void IrregularWaves::ReadEtaFromFile() {
    hydroc::debug::LogDebug(std::string("Reading eta file ") + params_.eta_file_path_ + ".");
    std::ifstream file(params_.eta_file_path_);
    if (!file) {
        throw std::runtime_error("Unable to open file at: " + params_.eta_file_path_ + ".");
    }

    std::string line;
    double      time;
    double      eta;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        char              delimiter;
        if (!(ss >> time >> delimiter >> eta) || delimiter != ':') {
            throw std::runtime_error("Could not parse line: " + line + ".");
        }
        time_data_.push_back(time);
        free_surface_elevation_sampled_.push_back(eta);
    }
    hydroc::debug::LogDebug("Finished reading eta file.");
}

Eigen::MatrixXd IrregularWaves::GetExcitationIRF(int b) const {
    return wave_info_[b].excitation_irf_matrix;
}

void IrregularWaves::AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data,
                               HydroData::SimulationParameters&           sim_data) {
    wave_info_   = irreg_h5_data;
    water_depth_ = sim_data.water_depth;
    g_           = sim_data.g;

    InitializeIRFVectors();
}

Eigen::Vector3d IrregularWaves::GetVelocity(const Eigen::Vector3d& position, double time, double elevation) {
    auto position_stretched =
        params_.wave_stretching_ ? GetWheelerStretchedPosition(position, elevation, water_depth_, mwl_) : position;
    return GetWaterVelocityIrregular(position_stretched, time, spectrum_frequencies_, spectral_densities_,
                                     spectral_widths_, wave_phases_, wavenumbers_, water_depth_, mwl_);
}

Eigen::Vector3d IrregularWaves::GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) {
    auto position_stretched =
        params_.wave_stretching_ ? GetWheelerStretchedPosition(position, elevation, water_depth_, mwl_) : position;
    return GetWaterAccelerationIrregular(position_stretched, time, spectrum_frequencies_, spectral_densities_,
                                         spectral_widths_, wave_phases_, wavenumbers_, water_depth_, mwl_);
}

double IrregularWaves::GetElevation(const Eigen::Vector3d& position, double time) {
    // Fallback for eta-file import mode (spectrum not generated).
    if (amplitudes_.size() == 0) {
        return GetEtaIrregular(position, time, spectrum_frequencies_, spectral_densities_,
                               spectral_widths_, wave_phases_, wavenumbers_);
    }

    // -------------------------------------------------------------------------
    // Free-surface elevation from linear superposition of wave components.
    //
    // The irregular sea state is represented as a sum of regular wave components:
    //
    //     η(x, t) = Σ A_i · cos(k_i·x − ω_i·t + φ_i)
    //
    // where:
    //     A_i   = wave amplitude for component i [m]
    //     k_i   = wavenumber [rad/m], from dispersion relation ω² = g·k·tanh(k·h)
    //     ω_i   = angular frequency [rad/s] = 2π·f_i
    //     φ_i   = random phase [rad], uniformly distributed in [0, 2π)
    //     x     = position along wave propagation direction [m]
    //     t     = time [s]
    //
    // Performance: Uses Eigen's vectorized operations (SIMD) for the summation.
    // The cos() is applied element-wise and the dot product sums the result.
    // -------------------------------------------------------------------------
    const double x = position.x();

    // Vectorized phase computation: phase_i = k_i*x - ω_i*t + φ_i
    const Eigen::ArrayXd phases = wavenumbers_.array() * x 
                                - angular_freqs_.array() * time 
                                + wave_phases_.array();

    // Vectorized elevation: η = Σ A_i * cos(phase_i)
    return (amplitudes_.array() * phases.cos()).sum();
}

Eigen::Vector2d IrregularWaves::GetElevationGradientXY(const Eigen::Vector3d& position, double time) const {
    // Fallback for eta-file import mode (spectrum not generated).
    if (amplitudes_.size() == 0) {
        double deta_dx = GetEtaGradientXIrregular(position, time, spectrum_frequencies_, spectral_densities_,
                                                   spectral_widths_, wave_phases_, wavenumbers_);
        return Eigen::Vector2d(deta_dx, 0.0);
    }

    // -------------------------------------------------------------------------
    // Free-surface slope (gradient) for visualization and normal computation.
    //
    // Taking the derivative of the elevation equation:
    //
    //     η(x, t)   = Σ A_i · cos(k_i·x − ω_i·t + φ_i)
    //     ∂η/∂x     = Σ −A_i · k_i · sin(k_i·x − ω_i·t + φ_i)
    //
    // Since waves propagate only in the +X direction, ∂η/∂y = 0.
    // The surface normal can be computed as: n = normalize(-∂η/∂x, -∂η/∂y, 1)
    //
    // Performance: Uses Eigen's vectorized operations (SIMD) for the summation.
    // -------------------------------------------------------------------------
    const double x = position.x();

    // Vectorized phase computation: phase_i = k_i*x - ω_i*t + φ_i
    const Eigen::ArrayXd phases = wavenumbers_.array() * x 
                                - angular_freqs_.array() * time 
                                + wave_phases_.array();

    // Vectorized gradient: ∂η/∂x = -Σ A_i * k_i * sin(phase_i)
    const double deta_dx = -(amplitudes_.array() * wavenumbers_.array() * phases.sin()).sum();

    return Eigen::Vector2d(deta_dx, 0.0);
}

double IrregularWaves::GetElevationForVisualization(const Eigen::Vector3d& position, 
                                                     double time, 
                                                     int max_components) const {
    // If no pre-computed amplitudes or max_components covers all, use full calculation.
    const Eigen::Index num_total = amplitudes_.size();
    if (num_total == 0) {
        return GetEtaIrregular(position, time, spectrum_frequencies_, spectral_densities_,
                               spectral_widths_, wave_phases_, wavenumbers_);
    }
    
    // Determine how many components to use.
    const Eigen::Index n = (max_components <= 0 || max_components >= num_total) 
                         ? num_total 
                         : static_cast<Eigen::Index>(max_components);
    
    // Use Eigen head() to get first n elements - still vectorized (SIMD).
    const double x = position.x();
    const Eigen::ArrayXd phases = wavenumbers_.head(n).array() * x 
                                - angular_freqs_.head(n).array() * time 
                                + wave_phases_.head(n).array();

    return (amplitudes_.head(n).array() * phases.cos()).sum();
}

Eigen::VectorXd IrregularWaves::GetForceAtTime(double t) {
    unsigned int total_dofs = num_bodies_ * 6;
    Eigen::VectorXd f(total_dofs);
    for (unsigned int i = 0; i < total_dofs; i++) {
        f[i] = 0.0;
    }

    for (unsigned int body = 0; body < num_bodies_; body++) {
        for (int dof = 0; dof < 6; ++dof) {
            double f_dof          = ExcitationConvolution(body, dof, t);
            unsigned int b_offset = body * 6;
            f[b_offset + dof]     = f_dof;
        }
    }

    return f;
}

Eigen::VectorXd IrregularWaves::SetSpectrumFrequencies(double start, double end, int num_points) {
    Eigen::VectorXd result(num_points);
    double step = (end - start) / (num_points - 1);

    for (int i = 0; i < num_points; ++i) {
        result[i] = start + i * step;
    }

    spectrum_frequencies_ = result;
    return result;
}

void IrregularWaves::CreateSpectrum() {
    int nf;
    if (params_.nfrequencies_ == 0) {
        double df = 1.0 / params_.simulation_duration_;
        nf        = std::ceil((params_.frequency_max_ - params_.frequency_min_) / df);
    } else {
        nf = static_cast<int>(params_.nfrequencies_);
    }
    spectrum_frequencies_ = Eigen::VectorXd::LinSpaced(nf, params_.frequency_min_, params_.frequency_max_);

    spectral_densities_ = JONSWAPSpectrumHz(spectrum_frequencies_, params_.wave_height_, params_.wave_period_,
                                            params_.peak_enhancement_factor_, params_.is_normalized_);

    spectral_widths_ = GetWidthArray(spectrum_frequencies_);

    wave_phases_ = Eigen::VectorXd(nf);
    std::mt19937 rng(params_.seed_);
    std::uniform_real_distribution<double> dist(0.0, 2 * M_PI);
    for (int i = 0; i < nf; ++i) {
        wave_phases_[i] = dist(rng);
    }

    auto omegas  = 2 * M_PI * spectrum_frequencies_;
    wavenumbers_ = ComputeWaveNumbers(omegas, water_depth_, g_);

    // Pre-compute amplitude and omega arrays for fast GetElevation().
    PrecomputeAmplitudes();
}

void IrregularWaves::PrecomputeAmplitudes() {
    // -------------------------------------------------------------------------
    // Pre-compute wave component amplitudes and angular frequencies.
    //
    // For each frequency component i:
    //     A_i     = sqrt(2 · S(f_i) · Δf_i)   [m]
    //     ω_i     = 2π · f_i                  [rad/s]
    //
    // where:
    //     S(f_i)  = spectral density at frequency f_i [m²/Hz]
    //     Δf_i    = frequency bin width [Hz]
    //
    // This pre-computation eliminates ~2N operations per GetElevation() call,
    // which is critical when evaluating thousands of grid points per frame.
    // -------------------------------------------------------------------------
    const Eigen::Index num_frequencies = spectrum_frequencies_.size();
    amplitudes_.resize(num_frequencies);
    angular_freqs_.resize(num_frequencies);

    for (Eigen::Index i = 0; i < num_frequencies; ++i) {
        // Wave amplitude from spectral density: A = sqrt(2 * S * df)
        amplitudes_[i] = std::sqrt(2.0 * spectral_densities_[i] * spectral_widths_[i]);
        
        // Angular frequency: omega = 2*pi*f
        angular_freqs_[i] = 2.0 * M_PI * spectrum_frequencies_[i];
    }
}

void IrregularWaves::CreateFreeSurfaceElevation() {
    double t_irf_min = 0.0;
    double t_irf_max = 0.0;
    for (size_t ii = 0; ii < ex_irf_time_sampled_.size(); ii++) {
        if (ex_irf_time_sampled_[ii][0] < t_irf_min) {
            t_irf_min = ex_irf_time_sampled_[ii][0];
        }
        if (ex_irf_time_sampled_[ii][0] > t_irf_max) {
            t_irf_max = ex_irf_time_sampled_[ii][0];
        }
        if (ex_irf_time_sampled_[ii][ex_irf_time_sampled_[ii].size() - 1] > t_irf_max) {
            t_irf_max = ex_irf_time_sampled_[ii][ex_irf_time_sampled_[ii].size() - 1];
        }
        if (ex_irf_time_sampled_[ii][ex_irf_time_sampled_[ii].size() - 1] < t_irf_min) {
            t_irf_min = ex_irf_time_sampled_[ii][ex_irf_time_sampled_[ii].size() - 1];
        }
    }
    auto duration      = params_.simulation_duration_ + 2 * (t_irf_max - t_irf_min);
    auto num_timesteps = static_cast<int>(std::ceil(duration / params_.simulation_dt_));
    auto time_array    = Eigen::VectorXd::LinSpaced(num_timesteps + 1, 0, num_timesteps * params_.simulation_dt_);

    free_surface_time_sampled_.resize(time_array.size());
    Eigen::VectorXd::Map(&free_surface_time_sampled_[0], time_array.size()) = time_array;
    for (size_t ii = 0; ii < free_surface_time_sampled_.size(); ii++) {
        free_surface_time_sampled_[ii] += -t_irf_max;
    }

    hydroc::debug::LogDebug(std::string("Precalculating free surface elevation from ") +
                            std::to_string(free_surface_time_sampled_.front()) + " to " +
                            std::to_string(free_surface_time_sampled_.back()) + ".");

    auto position = Eigen::Vector3d(0.0, 0.0, 0.0);
    free_surface_elevation_sampled_ =
        GetEtaIrregularTimeSeries(position, free_surface_time_sampled_, spectrum_frequencies_, spectral_densities_,
                                  spectral_widths_, wave_phases_, wavenumbers_);

    if (params_.ramp_duration_ > 0.0) {
        for (size_t i = 0; i < free_surface_time_sampled_.size(); ++i) {
            if (free_surface_time_sampled_[i] < params_.ramp_duration_) {
                if (free_surface_time_sampled_[i] <= 0.0) {
                    free_surface_elevation_sampled_[i] *= 0.0;
                } else {
                    free_surface_elevation_sampled_[i] *= free_surface_time_sampled_[i] / params_.ramp_duration_;
                }
            }
        }
    }

    hydroc::debug::LogDebug("Finished precalculating free surface elevation.");
}

void IrregularWaves::ResampleIRF(double dt) {
    for (unsigned int b = 0; b < num_bodies_; b++) {
        auto& time_array = ex_irf_time_sampled_[b];
        // Note: width_array is recalculated by CalculateWidthIRF() below
        auto& val_array  = ex_irf_sampled_[b];

        auto time_array_old = time_array;

        auto t0    = time_array_old[0];
        auto t1    = time_array_old[time_array_old.size() - 1];
        time_array = Eigen::VectorXd::LinSpaced(static_cast<int>(std::ceil((t1 - t0) / dt)), t0, t1);

        CalculateWidthIRF();

        assert(val_array.rows() == 6);
        Eigen::MatrixXd vals_new(6, time_array.size());

        Eigen::VectorXd t_old_scaled = Eigen::VectorXd::LinSpaced(time_array_old.size(), 0, 1);
        Eigen::VectorXd t_new_scaled = Eigen::VectorXd::LinSpaced(time_array.size(), 0, 1);

        Eigen::Spline<double, 6> spline =
            Eigen::SplineFitting<Eigen::Spline<double, 6>>::Interpolate(val_array, 3, t_old_scaled);
        for (int i = 0; i < time_array.rows(); i++) {
            vals_new.col(i) = spline(t_new_scaled[i]);
        }

        val_array = vals_new;
    }
}

void IrregularWaves::CalculateWidthIRF() {
    for (unsigned int b = 0; b < num_bodies_; b++) {
        auto& time_array  = ex_irf_time_sampled_[b];
        auto& width_array = ex_irf_width_sampled_[b];
        width_array       = GetWidthArray(time_array);
    }
}

// Return last index of vector element below value.
// - value: Input value
// - ticks: Array of ticks from which to find lower-bound index (assuming ascending order)
static size_t get_lower_index(double value, const std::vector<double>& ticks) {
    auto it = std::upper_bound(ticks.begin(), ticks.end(), value);
    // get nearest-below index
    size_t idx = it - ticks.begin() - 1;
    // remove one if equal to value
    if (ticks[idx] == value) {
        idx -= 1;
    }
    if (idx <= 0 || idx >= ticks.size() - 1) {
        throw std::runtime_error("Could not find index for value " + std::to_string(value) + " in array with bounds (" +
                                 std::to_string(ticks.front()) + ", " + std::to_string(ticks.back()) + ").");
    }
    // return index
    return idx;
}

double IrregularWaves::ExcitationConvolution(int body, int dof, double time) {
    double f_ex          = 0.0;
    auto&  irf_time_array  = ex_irf_time_sampled_[body];
    auto&  irf_val_mat     = ex_irf_sampled_[body];
    auto&  irf_width_array = ex_irf_width_sampled_[body];

    auto tmin = free_surface_time_sampled_.front();
    auto tmax = free_surface_time_sampled_.back();
    double t_tau_init = time - irf_time_array[0];
    int idx           = 0;
    if (t_tau_init <= tmin) {
        idx = 0;
    } else if (t_tau_init >= tmax) {
        idx = static_cast<int>(free_surface_time_sampled_.size()) - 2;
    } else {
        idx = static_cast<int>(get_lower_index(t_tau_init, free_surface_time_sampled_));
    }

    for (Eigen::Index j = 0; j < irf_time_array.size(); ++j) {
        double tau   = irf_time_array[j];
        double t_tau = time - tau;
        if (tmin <= t_tau && t_tau <= tmax) {
            while (free_surface_time_sampled_[idx] > t_tau) {
                idx -= 1;
            }

            auto t1 = free_surface_time_sampled_[idx];
            auto t2 = free_surface_time_sampled_[idx + 1];

            double eta_val;
            if (t_tau == t1) {
                eta_val = free_surface_elevation_sampled_[idx];
            } else if (t_tau == t2) {
                eta_val = free_surface_elevation_sampled_[idx + 1];
            } else if (t_tau > t1 && t_tau < t2) {
                auto eta1 = free_surface_elevation_sampled_[idx];
                auto eta2 = free_surface_elevation_sampled_[idx + 1];
                auto w1   = (t2 - t_tau) / (t2 - t1);
                auto w2   = 1.0 - w1;
                eta_val   = w1 * eta1 + w2 * eta2;
            } else {
                throw std::runtime_error("Excitation convolution: wrong tau value " + std::to_string(tau) +
                                         " not between " + std::to_string(t1) + " and " + std::to_string(t2) + ".");
            }

            f_ex += irf_val_mat(dof, j) * eta_val * irf_width_array[j];

        } else {
            // Time is outside precomputed range - skip this convolution term.
            // This can happen at the very end of simulation when time slightly exceeds
            // the simulation duration due to timestep discretization.
            continue;
        }
    }

    return f_ex;
}

