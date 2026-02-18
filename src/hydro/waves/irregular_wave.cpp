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
#include <hydroc/logging.h>

IrregularWaves::IrregularWaves(const IrregularWaveParams& params) : params_(params) {
    wave_stretching_ = params.wave_stretching;
}

void IrregularWaves::InitializeIRFVectors() {
    ex_irf_sampled_.resize(num_bodies_);
    ex_irf_time_sampled_.resize(num_bodies_);
    ex_irf_width_sampled_.resize(num_bodies_);

    for (unsigned int b = 0; b < num_bodies_; b++) {
        ex_irf_sampled_[b]      = GetExcitationIRF(b);
        ex_irf_time_sampled_[b] = wave_info_[b].excitation_irf_time;
    }
    CalculateWidthIRF();

    // Cache maximum IRF time span (used for ramp boundary check and eta padding).
    irf_time_max_ = 0.0;
    for (unsigned int b = 0; b < num_bodies_; b++) {
        if (ex_irf_time_sampled_[b].size() > 0) {
            double t_end = ex_irf_time_sampled_[b][ex_irf_time_sampled_[b].size() - 1];
            if (t_end > irf_time_max_) irf_time_max_ = t_end;
        }
    }

    if (!params_.eta_file_path.empty()) {
        std::vector<double> time_data = ReadEtaFromFile();

        if (time_data.empty() || free_surface_elevation_sampled_.empty()) {
            throw std::runtime_error("Failed to read eta data from file: " + params_.eta_file_path);
        }

        double eta_tmin = time_data.front();
        double required_tmin = -irf_time_max_;

        if (required_tmin < eta_tmin && time_data.size() > 1) {
            double dt = time_data[1] - time_data[0];
            int num_pad = static_cast<int>(std::ceil((eta_tmin - required_tmin) / dt));

            std::vector<double> extended_time(num_pad + time_data.size());
            std::vector<double> extended_eta(num_pad + free_surface_elevation_sampled_.size());

            for (int i = 0; i < num_pad; ++i) {
                extended_time[i] = eta_tmin - (num_pad - i) * dt;
                extended_eta[i] = 0.0;
            }
            for (size_t i = 0; i < time_data.size(); ++i) {
                extended_time[num_pad + i] = time_data[i];
                extended_eta[num_pad + i] = free_surface_elevation_sampled_[i];
            }

            free_surface_time_sampled_ = std::move(extended_time);
            free_surface_elevation_sampled_ = std::move(extended_eta);
        } else {
            free_surface_time_sampled_ = std::move(time_data);
        }

        use_eta_from_file_ = true;
        spectrum_created_ = false;
    } else if (params_.wave_height != 0.0 && params_.wave_period != 0.0) {
        CreateSpectrum();
        PrecomputeExcitationTransfer();
        use_eta_from_file_ = false;
        spectrum_created_ = true;
    }
}

std::vector<double> IrregularWaves::GetSpectrum() {
    if (!spectrum_created_) {
        throw std::runtime_error("Spectrum has not been created. Initialize with wave height and period to create spectrum.");
    }
    return std::vector<double>(spectral_densities_.data(),
                               spectral_densities_.data() + spectral_densities_.size());
}

const std::vector<double>& IrregularWaves::GetFreeSurfaceElevation() const {
    return free_surface_elevation_sampled_;
}

const std::vector<double>& IrregularWaves::GetFreeSurfaceTime() const {
    return free_surface_time_sampled_;
}

std::vector<double> IrregularWaves::GetFrequenciesHz() const {
    std::vector<double> out(spectrum_frequencies_.size());
    for (Eigen::Index i = 0; i < spectrum_frequencies_.size(); ++i) {
        out[i] = spectrum_frequencies_[i];
    }
    return out;
}

std::vector<double> IrregularWaves::ReadEtaFromFile() {
    hydroc::debug::LogDebug(std::string("Reading eta file ") + params_.eta_file_path + ".");
    std::ifstream file(params_.eta_file_path);
    if (!file) {
        throw std::runtime_error("Unable to open file at: " + params_.eta_file_path + ".");
    }

    std::vector<double> time_data;
    std::string line;
    double      time;
    double      eta;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        char              delimiter;
        if (!(ss >> time >> delimiter >> eta) || delimiter != ':') {
            throw std::runtime_error("Could not parse line: " + line + ".");
        }
        time_data.push_back(time);
        free_surface_elevation_sampled_.push_back(eta);
    }
    hydroc::debug::LogDebug("Finished reading eta file.");
    return time_data;
}

const Eigen::MatrixXd& IrregularWaves::GetExcitationIRF(int b) const {
    return wave_info_[b].excitation_irf_matrix;
}

void IrregularWaves::AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data,
                               const HydroData::SimulationParameters&     sim_data) {
    wave_info_   = irreg_h5_data;
    water_depth_ = sim_data.water_depth;
    g_           = sim_data.g;

    InitializeIRFVectors();
}

Eigen::Vector3d IrregularWaves::GetVelocity(const Eigen::Vector3d& position, double time, double elevation) const {
    auto position_stretched =
        params_.wave_stretching ? GetWheelerStretchedPosition(position, elevation, water_depth_, mwl_) : position;
    return GetWaterVelocityIrregular(position_stretched, time, spectrum_frequencies_, spectral_densities_,
                                     spectral_widths_, wave_phases_, wavenumbers_, water_depth_, mwl_);
}

Eigen::Vector3d IrregularWaves::GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) const {
    auto position_stretched =
        params_.wave_stretching ? GetWheelerStretchedPosition(position, elevation, water_depth_, mwl_) : position;
    return GetWaterAccelerationIrregular(position_stretched, time, spectrum_frequencies_, spectral_densities_,
                                         spectral_widths_, wave_phases_, wavenumbers_, water_depth_, mwl_);
}

double IrregularWaves::GetElevation(const Eigen::Vector3d& position, double time) const {
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

Eigen::VectorXd IrregularWaves::GetForceAtTime(double t) const {
    unsigned int total_dofs = num_bodies_ * kDofsPerBody;
    Eigen::VectorXd f = Eigen::VectorXd::Zero(total_dofs);

    if (use_eta_from_file_) {
        for (unsigned int body = 0; body < num_bodies_; body++) {
            unsigned int b_offset = body * kDofsPerBody;
            for (int dof = 0; dof < kDofsPerBody; ++dof) {
                f[b_offset + dof] = ExcitationConvolution(body, dof, t);
            }
        }
        return f;
    }

    bool ramp_active = (params_.ramp_duration > 0.0) && (t < params_.ramp_duration + irf_time_max_);

    // Phase terms shared by all bodies and DOFs at this timestep.
    const Eigen::ArrayXd theta = wave_phases_.array() - angular_freqs_.array() * t;
    const Eigen::ArrayXd cos_theta = theta.cos();
    const Eigen::ArrayXd sin_theta = theta.sin();

    for (unsigned int body = 0; body < num_bodies_; body++) {
        unsigned int b_offset = body * kDofsPerBody;

        if (ramp_active) {
            // Batch eta evaluation via matrix-vector product.
            // eta_j = cos_wt(j,:) · (A.*cos_theta) + sin_wt(j,:) · (A.*sin_theta)
            Eigen::VectorXd a_cos = (amplitudes_.array() * cos_theta).matrix();
            Eigen::VectorXd a_sin = (amplitudes_.array() * sin_theta).matrix();
            Eigen::VectorXd eta = irf_cos_wt_[body] * a_cos - irf_sin_wt_[body] * a_sin;

            auto& irf_time  = ex_irf_time_sampled_[body];
            auto& irf_val   = ex_irf_sampled_[body];
            auto& irf_width = ex_irf_width_sampled_[body];

            for (Eigen::Index j = 0; j < eta.size(); ++j) {
                double t_eval = t - irf_time[j];
                if (t_eval <= 0.0) {
                    eta[j] = 0.0;
                } else if (t_eval < params_.ramp_duration) {
                    eta[j] *= t_eval / params_.ramp_duration;
                }
            }

            Eigen::VectorXd eta_w = eta.array() * irf_width.array();
            for (int dof = 0; dof < kDofsPerBody; ++dof) {
                f[b_offset + dof] = irf_val.row(dof).dot(eta_w);
            }
        } else {
            // Fast path: frequency-domain transfer function.
            // C/S rows must be transposed to column arrays for element-wise ops.
            const auto& C = ex_transfer_cos_[body];
            const auto& S = ex_transfer_sin_[body];
            for (int dof = 0; dof < kDofsPerBody; ++dof) {
                Eigen::ArrayXd c_dof = C.row(dof).transpose().array();
                Eigen::ArrayXd s_dof = S.row(dof).transpose().array();
                f[b_offset + dof] = (amplitudes_.array() *
                    (c_dof * cos_theta - s_dof * sin_theta)).sum();
            }
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
    if (params_.nfrequencies <= 0) {
        throw std::invalid_argument("IrregularWaves::CreateSpectrum: nfrequencies must be > 0");
    }
    int nf = params_.nfrequencies;
    spectrum_frequencies_ = Eigen::VectorXd::LinSpaced(nf, params_.frequency_min, params_.frequency_max);

    spectral_densities_ = JONSWAPSpectrumHz(spectrum_frequencies_, params_.wave_height, params_.wave_period,
                                            params_.peak_enhancement_factor, params_.is_normalized);

    spectral_widths_ = GetWidthArray(spectrum_frequencies_);

    wave_phases_ = Eigen::VectorXd(nf);
    std::mt19937 rng(params_.seed);
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

void IrregularWaves::PrecomputeExcitationTransfer() {
    // -------------------------------------------------------------------------
    // Pre-compute matrices for fast excitation force evaluation.
    //
    // 1. Trig matrices: cos(ω_i τ_j) and sin(ω_i τ_j) for all (j, i).
    //    Used for batch eta evaluation during the ramp period via matrix-vector product.
    //
    // 2. Transfer function coefficients (derived from the trig matrices):
    //      C(dof,i) = Σ_j K(dof,j) Δτ_j cos(ω_i τ_j)
    //      S(dof,i) = Σ_j K(dof,j) Δτ_j sin(ω_i τ_j)
    //    These allow the convolution to be evaluated as:
    //      F_ex(dof,t) = Σ_i A_i [C(dof,i) cos(θ_i) - S(dof,i) sin(θ_i)]
    //    reducing the per-timestep cost from O(N_irf × N_freq) to O(N_freq).
    // -------------------------------------------------------------------------
    const Eigen::Index nf = angular_freqs_.size();

    irf_cos_wt_.resize(num_bodies_);
    irf_sin_wt_.resize(num_bodies_);
    ex_transfer_cos_.resize(num_bodies_);
    ex_transfer_sin_.resize(num_bodies_);

    for (unsigned int b = 0; b < num_bodies_; ++b) {
        auto& irf_time  = ex_irf_time_sampled_[b];
        auto& irf_val   = ex_irf_sampled_[b];
        auto& irf_width = ex_irf_width_sampled_[b];
        const Eigen::Index n_irf = irf_time.size();

        irf_cos_wt_[b].resize(n_irf, nf);
        irf_sin_wt_[b].resize(n_irf, nf);

        for (Eigen::Index j = 0; j < n_irf; ++j) {
            for (Eigen::Index i = 0; i < nf; ++i) {
                double omega_tau = angular_freqs_[i] * irf_time[j];
                irf_cos_wt_[b](j, i) = std::cos(omega_tau);
                irf_sin_wt_[b](j, i) = std::sin(omega_tau);
            }
        }

        // Transfer function: C = (K.*dτ)^T * cos_matrix, S = (K.*dτ)^T * sin_matrix
        ex_transfer_cos_[b].resize(kDofsPerBody, nf);
        ex_transfer_sin_[b].resize(kDofsPerBody, nf);
        for (int dof = 0; dof < kDofsPerBody; ++dof) {
            Eigen::VectorXd kw = irf_val.row(dof).transpose().array() * irf_width.array();
            ex_transfer_cos_[b].row(dof) = kw.transpose() * irf_cos_wt_[b];
            ex_transfer_sin_[b].row(dof) = kw.transpose() * irf_sin_wt_[b];
        }

    }
}

std::pair<std::vector<double>, std::vector<double>>
IrregularWaves::ComputeElevationTimeSeries(double t_start, double t_end, double dt) const {
    int n = static_cast<int>(std::ceil((t_end - t_start) / dt)) + 1;
    std::vector<double> times(n);
    std::vector<double> elevations(n);
    const Eigen::Vector3d origin(0.0, 0.0, 0.0);

    for (int i = 0; i < n; ++i) {
        double t = t_start + i * dt;
        times[i] = t;
        elevations[i] = GetElevation(origin, t);
    }
    return {times, elevations};
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

    auto dist = static_cast<ptrdiff_t>(it - ticks.begin());
    if (dist <= 0) {
        throw std::runtime_error("Could not find index for value " + std::to_string(value) +
                                 " in array with bounds (" + std::to_string(ticks.front()) +
                                 ", " + std::to_string(ticks.back()) + ").");
    }
    auto idx = static_cast<size_t>(dist - 1);

    if (idx > 0 && ticks[idx] == value) {
        --idx;
    }
    if (idx == 0 || idx >= ticks.size() - 1) {
        throw std::runtime_error("Could not find index for value " + std::to_string(value) +
                                 " in array with bounds (" + std::to_string(ticks.front()) +
                                 ", " + std::to_string(ticks.back()) + ").");
    }
    return idx;
}

double IrregularWaves::ExcitationConvolution(int body, int dof, double time) const {
    // Eta file import path only — spectrum-based forces are computed in GetForceAtTime().
    double f_ex = 0.0;
    auto& irf_time_array  = ex_irf_time_sampled_[body];
    auto& irf_val_mat     = ex_irf_sampled_[body];
    auto& irf_width_array = ex_irf_width_sampled_[body];

    auto tmin = free_surface_time_sampled_.front();
    auto tmax = free_surface_time_sampled_.back();
    double t_tau_init = time - irf_time_array[0];
    int idx = 0;
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
                throw std::runtime_error(
                    "Excitation convolution: wrong tau value " + std::to_string(tau) +
                    " not between " + std::to_string(t1) + " and " + std::to_string(t2) + ".");
            }

            f_ex += irf_val_mat(dof, j) * eta_val * irf_width_array[j];
        }
    }
    return f_ex;
}

