/*********************************************************************
 * @file  regular_wave.cpp
 * @brief Regular wave model implementation.
 *********************************************************************/

#include <hydroc/waves/regular_wave.h>
#include "wave_utilities.h"

#include <cmath>

RegularWave::RegularWave()
    : num_bodies_(1),
      wavenumber_(0.0) {
}

RegularWave::RegularWave(unsigned int num_b)
    : num_bodies_(num_b),
      wavenumber_(0.0) {
}

RegularWave::RegularWave(const RegularWaveParams& params) {
    num_bodies_ = params.num_bodies_;
    regular_wave_amplitude_ = params.regular_wave_amplitude_;
    regular_wave_omega_     = params.regular_wave_omega_;
    regular_wave_phase_     = params.regular_wave_phase_;
    wave_stretching_        = params.wave_stretching_;
}

void RegularWave::Initialize() {
    wavenumber_ = ComputeWaveNumber(regular_wave_omega_, water_depth_, g_);
}

void RegularWave::AddH5Data(std::vector<HydroData::RegularWaveInfo>& reg_h5_data, HydroData::SimulationParameters& sim_data) {
    wave_info_   = reg_h5_data;
    water_depth_ = sim_data.water_depth;
    g_           = sim_data.g;

    int total_dofs = 6 * num_bodies_;
    excitation_force_mag_.resize(total_dofs);
    excitation_force_phase_.resize(total_dofs);
    force_.resize(total_dofs);

    double wave_omega_delta = GetOmegaDelta();
    double freq_index_des   = (regular_wave_omega_ / wave_omega_delta) - 1;
    for (unsigned int b = 0; b < num_bodies_; b++) {
        for (int rowEx = 0; rowEx < 6; rowEx++) {
            unsigned int body_offset = 6 * b;
            excitation_force_mag_[body_offset + rowEx]   = GetExcitationMagInterp(b, rowEx, 0, freq_index_des);
            excitation_force_phase_[body_offset + rowEx] = GetExcitationPhaseInterp(b, rowEx, 0, freq_index_des);
        }
    }
}

Eigen::Vector3d RegularWave::GetVelocity(const Eigen::Vector3d& position, double time, double elevation) {
    auto position_stretched =
        wave_stretching_ ? GetWheelerStretchedPosition(position, elevation, water_depth_, mwl_) : position;
    return GetWaterVelocity(position_stretched, time, regular_wave_omega_, regular_wave_amplitude_, regular_wave_phase_,
                            wavenumber_, water_depth_, mwl_);
}

Eigen::Vector3d RegularWave::GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) {
    auto position_stretched =
        wave_stretching_ ? GetWheelerStretchedPosition(position, elevation, water_depth_, mwl_) : position;
    return GetWaterAcceleration(position_stretched, time, regular_wave_omega_, regular_wave_amplitude_,
                                regular_wave_phase_, wavenumber_, water_depth_, mwl_);
}
  
double RegularWave::GetElevation(const Eigen::Vector3d& position, double time) {
    return GetEta(position, time, regular_wave_omega_, regular_wave_amplitude_, regular_wave_phase_, wavenumber_);
}

Eigen::Vector2d RegularWave::GetElevationGradientXY(const Eigen::Vector3d& position, double time) const {
    // Regular waves propagate in +X direction only; ∂η/∂y = 0.
    double deta_dx = GetEtaGradientX(position, time, regular_wave_omega_, regular_wave_amplitude_,
                                      regular_wave_phase_, wavenumber_);
    return Eigen::Vector2d(deta_dx, 0.0);
}

Eigen::VectorXd RegularWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies_ * 6;
    Eigen::VectorXd f(dof);
    for (unsigned int b = 0; b < num_bodies_; b++) {
        unsigned int body_offset = 6 * b;
        for (int rowEx = 0; rowEx < 6; rowEx++) {
            f[body_offset + rowEx] = excitation_force_mag_[body_offset + rowEx] * regular_wave_amplitude_ *
                                     std::cos(regular_wave_omega_ * t + excitation_force_phase_[body_offset + rowEx]);
        }
    }
    return f;
}

double RegularWave::GetOmegaDelta() const {
    double omega_max = wave_info_[0].freq_list[wave_info_[0].freq_list.size() - 1];
    double num_freqs = wave_info_[0].freq_list.size();
    return omega_max / num_freqs;
}

double RegularWave::GetExcitationMagInterp(int b, int i, int j, double freq_index_des) const {
    double freq_interp_val    = freq_index_des - std::floor(freq_index_des);
    double excitationMagFloor = wave_info_[b].excitation_mag_matrix(i, j, (int)std::floor(freq_index_des));
    double excitationMagCeil  = wave_info_[b].excitation_mag_matrix(i, j, (int)std::floor(freq_index_des) + 1);
    double excitationMag      = (freq_interp_val * (excitationMagCeil - excitationMagFloor)) + excitationMagFloor;

    return excitationMag;
}

double RegularWave::GetExcitationPhaseInterp(int b, int i, int j, double freq_index_des) const {
    double freq_interp_val      = freq_index_des - std::floor(freq_index_des);
    double excitationPhaseFloor = wave_info_[b].excitation_phase_matrix(i, j, (int)std::floor(freq_index_des));
    double excitationPhaseCeil  = wave_info_[b].excitation_phase_matrix(i, j, (int)std::floor(freq_index_des) + 1);
    double excitationPhase =
        (freq_interp_val * (excitationPhaseCeil - excitationPhaseFloor)) + excitationPhaseFloor;

    return excitationPhase;
}

