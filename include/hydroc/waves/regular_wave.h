/*********************************************************************
 * @file  regular_wave.h
 * @brief Regular wave model declarations.
 *********************************************************************/

#ifndef HYDROC_WAVES_REGULAR_WAVE_H
#define HYDROC_WAVES_REGULAR_WAVE_H

#include <hydroc/waves/wave_base.h>
#include <hydroc/io/h5_reader.h>

#include <vector>

// Forward declaration for internal utilities (not exposed in public API)
// Wave utilities are internal implementation details

struct RegularWaveParams {
    unsigned int num_bodies_;
    double regular_wave_amplitude_;
    double regular_wave_omega_;
    double regular_wave_phase_ = 0.0;
    bool wave_stretching_      = true;
};

class RegularWave : public WaveBase {
  public:
    RegularWave();
    explicit RegularWave(unsigned int num_b);
    explicit RegularWave(const RegularWaveParams& params);

    void Initialize() override;
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode_; }

    //// RADU - eliminate these and use a RegularWaveParams struct (consistent with IrregularWaves)
    double regular_wave_amplitude_;
    double regular_wave_omega_;
    double regular_wave_phase_ = 0.0;

    void AddH5Data(std::vector<HydroData::RegularWaveInfo>& reg_h5_data, HydroData::SimulationParameters& sim_data);

    double GetElevation(const Eigen::Vector3d& position, double time) override;
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time, double elevation) override;
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) override;

  private:
    unsigned int num_bodies_ = 0;
    const WaveMode mode_ = WaveMode::regular;
    std::vector<HydroData::RegularWaveInfo> wave_info_;
    Eigen::VectorXd excitation_force_mag_;
    Eigen::VectorXd excitation_force_phase_;
    Eigen::VectorXd force_;
    double wavenumber_ = 0.0;

    double GetOmegaDelta() const;
    double GetExcitationMagInterp(int b, int i, int j, double freq_index_des) const;
    double GetExcitationPhaseInterp(int b, int i, int j, double freq_index_des) const;
};

#endif  // HYDROC_WAVES_REGULAR_WAVE_H

