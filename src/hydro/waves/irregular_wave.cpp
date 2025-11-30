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

#include <hydroc/helper.h>
#include <hydroc/logging.h>

IrregularWaves::IrregularWaves(const IrregularWaveParams& params) : params_(params) {}

void IrregularWaves::InitializeIRFVectors() {
    ex_irf_sampled_.resize(params_.num_bodies_);
    ex_irf_time_sampled_.resize(params_.num_bodies_);
    ex_irf_width_sampled_.resize(params_.num_bodies_);

    for (unsigned int b = 0; b < params_.num_bodies_; b++) {
        ex_irf_sampled_[b]      = GetExcitationIRF(b);
        ex_irf_time_sampled_[b] = wave_info_[b].excitation_irf_time;
        CalculateWidthIRF();
    }

    if (params_.simulation_dt_ > 0.0) {
        ResampleIRF(params_.simulation_dt_);
    }

    if (!params_.eta_file_path_.empty()) {
        ReadEtaFromFile();
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
    for (int i = 0; i < spectrum_frequencies_.size(); ++i) {
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

Eigen::Vector3d IrregularWaves::GetVelocity(const Eigen::Vector3d& position, double time) {
    auto position_stretched = position;
    if (params_.wave_stretching_) {
        auto eta = GetEtaIrregular(position, time, spectrum_frequencies_, spectral_densities_, spectral_widths_, wave_phases_, wavenumbers_);
        auto z_pos = position.z() - mwl_;
        position_stretched[2] = water_depth_ * (z_pos - eta) / (water_depth_ + eta);
    }

    return GetWaterVelocityIrregular(position_stretched,
                                     time,
                                     spectrum_frequencies_,
                                     spectral_densities_,
                                     spectral_widths_,
                                     wave_phases_,
                                     wavenumbers_,
                                     water_depth_,
                                     mwl_);
}

Eigen::Vector3d IrregularWaves::GetAcceleration(const Eigen::Vector3d& position, double time) {
    auto position_stretched = position;
    if (params_.wave_stretching_) {
        auto eta = GetEtaIrregular(position, time, spectrum_frequencies_, spectral_densities_, spectral_widths_, wave_phases_, wavenumbers_);
        auto z_pos = position.z() - mwl_;
        position_stretched[2] = water_depth_ * (z_pos - eta) / (water_depth_ + eta);
    }

    return GetWaterAccelerationIrregular(position_stretched,
                                         time,
                                         spectrum_frequencies_,
                                         spectral_densities_,
                                         spectral_widths_,
                                         wave_phases_,
                                         wavenumbers_,
                                         water_depth_,
                                         mwl_);
}

double IrregularWaves::GetElevation(const Eigen::Vector3d& position, double time) {
    return GetEtaIrregular(position, time, spectrum_frequencies_, spectral_densities_, spectral_widths_, wave_phases_, wavenumbers_);
}

Eigen::VectorXd IrregularWaves::GetForceAtTime(double t) {
    unsigned int total_dofs = params_.num_bodies_ * 6;
    Eigen::VectorXd f(total_dofs);
    for (int i = 0; i < total_dofs; i++) {
        f[i] = 0.0;
    }

    for (int body = 0; body < params_.num_bodies_; body++) {
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

void IrregularWaves::SetUpWaveMesh(std::string filename) {
    mesh_file_name_   = filename;
    int num_timesteps = static_cast<int>(std::ceil(params_.simulation_duration_ / params_.simulation_dt_));
    Eigen::VectorXd time_index = Eigen::VectorXd::LinSpaced(num_timesteps + 1, 0, num_timesteps * params_.simulation_dt_);
    Eigen::VectorXd eta_vec    = Eigen::Map<const Eigen::VectorXd>(free_surface_elevation_sampled_.data(),
                                                                   static_cast<Eigen::Index>(free_surface_elevation_sampled_.size()));

    std::vector<std::array<double, 3>> free_surface_3d_pts = CreateFreeSurface3DPts(eta_vec, time_index);
    std::vector<std::array<size_t, 3>> free_surface_triangles =
        CreateFreeSurfaceTriangles(static_cast<size_t>(time_index.size()));

    WriteFreeSurfaceMeshObj(free_surface_3d_pts, free_surface_triangles, mesh_file_name_);
}

std::string IrregularWaves::GetMeshFile() {
    return mesh_file_name_;
}

Eigen::Vector3<double> IrregularWaves::GetWaveMeshVelocity() {
    return Eigen::Vector3d(1.0, 0, 0);
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
}

void IrregularWaves::CreateFreeSurfaceElevation() {
    double t_irf_min = 0.0;
    double t_irf_max = 0.0;
    for (auto ii = 0; ii < ex_irf_time_sampled_.size(); ii++) {
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
    for (int ii = 0; ii < free_surface_time_sampled_.size(); ii++) {
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
    for (unsigned int b = 0; b < params_.num_bodies_; b++) {
        auto& time_array  = ex_irf_time_sampled_[b];
        auto& width_array = ex_irf_width_sampled_[b];
        auto& val_array   = ex_irf_sampled_[b];

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
    for (unsigned int b = 0; b < params_.num_bodies_; b++) {
        auto& time_array  = ex_irf_time_sampled_[b];
        auto& width_array = ex_irf_width_sampled_[b];
        width_array       = GetWidthArray(time_array);
    }
}

double IrregularWaves::ExcitationConvolution(int body, int dof, double time) {
    double f_ex          = 0.0;
    auto&  irf_time_array  = ex_irf_time_sampled_[body];
    auto&  irf_val_mat     = ex_irf_sampled_[body];
    auto&  irf_width_array = ex_irf_width_sampled_[body];

    auto tmin = free_surface_time_sampled_.front();
    auto tmax = free_surface_time_sampled_.back();
    double t_tau = time - irf_time_array[0];
    int idx      = 0;
    if (t_tau <= tmin) {
        idx = 0;
    } else if (t_tau >= tmax) {
        idx = static_cast<int>(free_surface_time_sampled_.size()) - 2;
    } else {
        idx = get_lower_index(t_tau, free_surface_time_sampled_);
    }

    for (size_t j = 0; j < irf_time_array.size(); ++j) {
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
            throw std::runtime_error(
                "Excitation convolution: trying to find free surface elevation at a time out of bounds from the "
                "precomputed free surface elevation (" +
                std::to_string(t_tau) + "not in [" + std::to_string(tmin) + ", " + std::to_string(tmax) +
                "]). Excitation force ignored at this time step.");
        }
    }

    return f_ex;
}

