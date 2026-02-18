/*********************************************************************
 * @file  irregular_wave.h
 * @brief Irregular wave model declarations.
 *********************************************************************/

#ifndef HYDROC_WAVES_IRREGULAR_WAVE_H
#define HYDROC_WAVES_IRREGULAR_WAVE_H

#include <hydroc/waves/wave_base.h>
#include <hydroc/io/h5_reader.h>

#include <string>
#include <vector>

struct IrregularWaveParams {
    double simulation_dt_ = 0.0;
    double simulation_duration_ = 0.0;
    double ramp_duration_ = 0.0;
    std::string eta_file_path_;
    double wave_height_             = 0.0;
    double wave_period_             = 0.0;
    double frequency_min_           = 0.001;
    double frequency_max_           = 1.0;
    double nfrequencies_            = 0;
    double peak_enhancement_factor_ = 1.0;
    bool is_normalized_             = false;
    int seed_                       = 1;
    bool wave_stretching_           = true;
};

class IrregularWaves : public WaveBase {
  public:
    explicit IrregularWaves(const IrregularWaveParams& params);
    void Initialize() override {}

    void CreateSpectrum();
    std::vector<double> GetSpectrum();
    std::vector<double> GetFreeSurfaceElevation();
    std::vector<double> GetFreeSurfaceTime() const;
    std::vector<double> GetFrequenciesHz() const;

    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode_; }

    Eigen::VectorXd SetSpectrumFrequencies(double start, double end, int num_steps);

    void AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data, HydroData::SimulationParameters& sim_data);

    double GetElevation(const Eigen::Vector3d& position, double time) override;
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time, double elevation) override;
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) override;

    /// Return the surface slope (∂η/∂x, ∂η/∂y) at a given position and time.
    /// Used for computing surface normals in visualization.
    Eigen::Vector2d GetElevationGradientXY(const Eigen::Vector3d& position, double time) const;

    /// Compute elevation using only the first max_components wave components.
    /// Used for faster visualization rendering when full accuracy is not required.
    /// @param max_components Maximum number of frequency components to use (0 or negative = use all)
    double GetElevationForVisualization(const Eigen::Vector3d& position, double time, int max_components) const;

  private:
    IrregularWaveParams params_;
    std::vector<double> spectrum_;
    std::vector<double> time_data_;
    std::vector<double> free_surface_elevation_sampled_;
    std::vector<double> free_surface_time_sampled_;
    bool spectrumCreated_ = false;

    const WaveMode mode_ = WaveMode::irregular;
    std::vector<HydroData::IrregularWaveInfo> wave_info_;
    std::vector<Eigen::MatrixXd> ex_irf_sampled_;
    std::vector<Eigen::VectorXd> ex_irf_time_sampled_;
    std::vector<Eigen::VectorXd> ex_irf_width_sampled_;
    Eigen::VectorXd spectrum_frequencies_;
    Eigen::VectorXd spectral_densities_;
    Eigen::VectorXd spectral_widths_;
    Eigen::VectorXd wavenumbers_;
    Eigen::VectorXd wave_phases_;

    // -------------------------------------------------------------------------
    // Pre-computed arrays for fast elevation calculation.
    // These are computed once from the spectrum and reused for every GetElevation() call.
    // This optimization is critical for real-time visualization where GetElevation()
    // may be called thousands of times per frame (once per water surface grid point).
    // -------------------------------------------------------------------------
    
    /// Wave amplitude for each frequency component [m].
    /// Computed as: A_i = sqrt(2 * S(f_i) * df_i), where S is the spectral density.
    Eigen::VectorXd amplitudes_;
    
    /// Angular frequency for each component [rad/s].
    /// Computed as: omega_i = 2 * pi * f_i
    Eigen::VectorXd angular_freqs_;

    void InitializeIRFVectors();
    
    /// Pre-compute amplitudes_ and angular_freqs_ arrays from the spectrum.
    /// Called automatically after CreateSpectrum(). Must be called before GetElevation().
    void PrecomputeAmplitudes();
    void ReadEtaFromFile();
    void CreateFreeSurfaceElevation();

    Eigen::MatrixXd GetExcitationIRF(int b) const;
    void ResampleIRF(double dt);
    void CalculateWidthIRF();
    double ExcitationConvolution(int body, int dof, double time);
};

#endif  // HYDROC_WAVES_IRREGULAR_WAVE_H

