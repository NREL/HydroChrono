/*********************************************************************
 * @file  regular_wave.h
 * @brief Regular wave model declarations.
 *********************************************************************/

#ifndef HYDRO_WAVES_REGULAR_WAVE_H
#define HYDRO_WAVES_REGULAR_WAVE_H

#include "wave_base.h"
#include "wave_utilities.h"

#include <hydroc/io/h5_reader.h>

#include <vector>

class RegularWave : public WaveBase {
  public:
    RegularWave();
    explicit RegularWave(unsigned int num_b);

    void Initialize() override;
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode_; }

    double regular_wave_amplitude_;
    double regular_wave_omega_;
    double regular_wave_phase_ = 0.0;

    void AddH5Data(std::vector<HydroData::RegularWaveInfo>& reg_h5_data, HydroData::SimulationParameters& sim_data);

    double GetElevation(const Eigen::Vector3d& position, double time) override;
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time) override;
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time) override;

  private:
    unsigned int num_bodies_;
    const WaveMode mode_ = WaveMode::regular;
    std::vector<HydroData::RegularWaveInfo> wave_info_;
    Eigen::VectorXd excitation_force_mag_;
    Eigen::VectorXd excitation_force_phase_;
    Eigen::VectorXd force_;
    double wavenumber_;

    double GetOmegaDelta() const;
    double GetExcitationMagInterp(int b, int i, int j, double freq_index_des) const;
    double GetExcitationPhaseInterp(int b, int i, int j, double freq_index_des) const;
};

#endif  // HYDRO_WAVES_REGULAR_WAVE_H

