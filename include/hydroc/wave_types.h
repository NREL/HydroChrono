#pragma once
#include <Eigen/Dense>

enum class WaveMode {
    /// @brief No waves
    noWaveCIC = 0,
    /// @brief Regular waves
    regular = 1,
    /// @brief Irregular waves
    irregular = 2
};

// pure virtual (interface) class for wave modes (regular, irregular, etc)
// use only Eigen3 types
class WaveBase {
  public:
    virtual void Initialize() = 0;
    virtual Eigen::VectorXd GetForceAtTime(double t) = 0;
    virtual WaveMode GetWaveMode() = 0;
};

// class to intstantiate WaveBase for no waves
class NoWave : public WaveBase {
  public:
    NoWave() { num_bodies = 1; }
    NoWave(unsigned int num_b) { num_bodies = num_b; }
    void Initialize() override { }
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode; }

  private:
    unsigned int num_bodies;
    const WaveMode mode = WaveMode::noWaveCIC;
};

// class to instantiate WaveBase for regular waves
class RegularWave : public WaveBase {
  public:
    RegularWave();
    RegularWave(unsigned int num_b);
    void Initialize() override;
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode; }

  private:
    unsigned int num_bodies;
    const WaveMode mode = WaveMode::regular;
};

// class to instantiate WaveBase for irregular waves
class IrregularWave : public WaveBase {
  public:
    IrregularWave();
    IrregularWave(unsigned int num_b);
    void Initialize() override;  // call any set up functions from here
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode; }
    // add more helper functions for calculations here:

  private:
    unsigned int num_bodies;
    const WaveMode mode = WaveMode::irregular;
};

// =============================================================================
class HydroInputs {
    WaveMode mode;
    HydroInputs();
    void UpdateNumTimesteps();
    void UpdateRampTimesteps();
    void CreateSpectrum();
    void CreateFreeSurfaceElevation();
    std::vector<double> spectrum_frequencies;
    std::vector<double> spectral_densities;
    std::vector<double> eta;
    double ramp_duration;
    double freq_index_des;
    double wave_height;
    double wave_period;
    double simulation_duration;
    double simulation_dt;
    int num_timesteps;
    int ramp_timesteps;
    std::vector<double> ramp;
    double regular_wave_amplitude;
    double regular_wave_omega;
    double wave_omega_delta;
    std::vector<double> excitation_force_mag;
    std::vector<double> excitation_force_phase;
    HydroInputs(HydroInputs& old) = default;
    HydroInputs& operator=(const HydroInputs& rhs) = default;
};
