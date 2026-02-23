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

/// Parameters for constructing a RegularWave.
///
/// Specify the wave frequency as either omega (rad/s) or period (s), but not
/// both.  If both are non-zero, the constructor throws std::invalid_argument.
struct RegularWaveParams {
    double regular_wave_amplitude = 0.0;
    double regular_wave_omega     = 0.0;
    double regular_wave_period    = 0.0;
    double regular_wave_phase     = 0.0;
    bool wave_stretching          = true;
};

class RegularWave : public WaveBase {
  public:
    RegularWave();
    explicit RegularWave(const RegularWaveParams& params);

    void Initialize() override;
    Eigen::VectorXd GetForceAtTime(double t) const override;
    WaveMode GetWaveMode() const override { return mode_; }

    /// @name Mutators
    /// @{
    void SetAmplitude(double amplitude);
    void SetOmega(double omega);
    void SetPeriod(double period);
    void SetPhase(double phase);
    /// @}

    /// @name Accessors
    /// @{
    double GetAmplitude() const { return regular_wave_amplitude_; }
    double GetOmega() const { return regular_wave_omega_; }
    double GetPeriod() const;
    double GetPhase() const { return regular_wave_phase_; }
    /// @}

    void AddH5Data(std::vector<HydroData::RegularWaveInfo>& reg_h5_data, const HydroData::SimulationParameters& sim_data);

    double GetElevation(const Eigen::Vector3d& position, double time) const override;
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time, double elevation) const override;
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) const override;

    /// Return the surface slope (∂η/∂x, ∂η/∂y) at a given position and time.
    /// Used for computing surface normals in visualization.
    Eigen::Vector2d GetElevationGradientXY(const Eigen::Vector3d& position, double time) const;

  private:
    static constexpr WaveMode mode_ = WaveMode::regular;

    double regular_wave_amplitude_ = 0.0;
    double regular_wave_omega_     = 0.0;
    double regular_wave_phase_     = 0.0;

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

