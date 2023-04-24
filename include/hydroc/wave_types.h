#pragma once
#include <hydroc/h5fileinfo.h>
#include <Eigen/Dense>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// todo move this helper function somewhere else?
Eigen::VectorXd PiersonMoskowitzSpectrumHz(Eigen::VectorXd& f, double Hs, double Tp);
Eigen::VectorXd FreeSurfaceElevation(const Eigen::VectorXd& freqs_hz,
                                     const Eigen::VectorXd& spectral_densities,
                                     const Eigen::VectorXd& time_index,
                                     int seed = 1);

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
    virtual void Initialize()                        = 0;
    virtual Eigen::VectorXd GetForceAtTime(double t) = 0;
    virtual WaveMode GetWaveMode()                   = 0;
};

// class to intstantiate WaveBase for no waves
class NoWave : public WaveBase {
  public:
    NoWave() { num_bodies = 1; }
    NoWave(unsigned int num_b) { num_bodies = num_b; }
    void Initialize() override {}
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

    // user input variables
    double regular_wave_amplitude;
    double regular_wave_omega;
    // double freq_index_des;
    // double wave_omega_delta;
    // Eigen::VectorXd excitation_force_mag;
    // Eigen::VectorXd excitation_force_phase;

    void AddH5Data(std::vector<HydroData::RegularWaveInfo>& reg_h5_data);

  private:
    unsigned int num_bodies;  // TODO is this needed?
    const WaveMode mode = WaveMode::regular;
    std::vector<HydroData::RegularWaveInfo> info;
    Eigen::VectorXd excitation_force_mag;
    Eigen::VectorXd excitation_force_phase;
    Eigen::VectorXd force;
    double GetOmegaDelta() const;
    double GetExcitationMagInterp(int b, int i, int j, double freq_index_des) const;
    double GetExcitationPhaseInterp(int b, int i, int j, double freq_index_des) const;
};

// class to instantiate WaveBase for irregular waves
class IrregularWave : public WaveBase {
  public:
    IrregularWave();
    IrregularWave(unsigned int num_b);
    void Initialize() override;  // call any set up functions from here
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode; }
    Eigen::VectorXd SetSpectrumFrequencies(double start, double end, int num_steps);
    void SetUpWaveMesh(std::string filename = "fse_mesh.obj");
    std::string GetMeshFile();
    Eigen::Vector3<double> GetWaveMeshVelocity();
    // add more helper functions for calculations here:

    double wave_height;
    double wave_period;
    double simulation_duration;
    double simulation_dt;
    double ramp_duration;
    Eigen::VectorXd eta;

    void AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data);

  private:
    unsigned int num_bodies;
    const WaveMode mode = WaveMode::irregular;
    std::vector<HydroData::IrregularWaveInfo> info;
    std::vector<Eigen::MatrixXd> ex_irf_resampled;
    std::vector<Eigen::VectorXd> ex_irf_time_resampled;
    Eigen::VectorXd spectrum_frequencies;
    Eigen::VectorXd spectral_densities;
    std::string mesh_file_name;

    Eigen::MatrixXd GetExcitationIRF(int b) const;
    Eigen::VectorXd ResampleTime(const Eigen::VectorXd& t_old, const double dt_new);
    Eigen::MatrixXd ResampleVals(const Eigen::VectorXd& t_old, Eigen::MatrixXd& vals_old, const Eigen::VectorXd& t_new);
    double ExcitationConvolution(int body, int dof, double time);

    void CreateSpectrum();
    void IrregularWave::CreateFreeSurfaceElevation();
    friend Eigen::VectorXd PiersonMoskowitzSpectrumHz(Eigen::VectorXd& f, double Hs, double Tp);
};

void WriteFreeSurfaceMeshObj(const std::vector<std::array<double, 3>>& points,
                             const std::vector<std::array<size_t, 3>>& triangles,
                             const std::string& file_name);
std::vector<std::array<size_t, 3>> CreateFreeSurfaceTriangles(size_t eta_size);
std::vector<std::array<double, 3>> CreateFreeSurface3DPts(const Eigen::VectorXd& eta, const Eigen::VectorXd& t_vec);
    // =============================================================================
// class HydroInputs {
//  public:
//    WaveMode mode;
//    HydroInputs();
//    void UpdateNumTimesteps();
//    void UpdateRampTimesteps();
//    void CreateSpectrum();
//    void CreateFreeSurfaceElevation();
//    std::vector<double> spectrum_frequencies;
//    std::vector<double> spectral_densities;
//    std::vector<double> eta;
//    double ramp_duration;
//    double freq_index_des;
//    double wave_height;
//    double wave_period;
//    double simulation_duration;
//    double simulation_dt;
//    int num_timesteps;
//    int ramp_timesteps;
//    std::vector<double> ramp;
//    double regular_wave_amplitude;
//    double regular_wave_omega;
//    double wave_omega_delta;
//    std::vector<double> excitation_force_mag;
//    std::vector<double> excitation_force_phase;
//    HydroInputs(HydroInputs& old) = default;
//    HydroInputs& operator=(const HydroInputs& rhs) = default;
//};
