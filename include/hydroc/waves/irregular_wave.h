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
    double wave_height             = 0.0;
    double wave_period             = 0.0;
    double frequency_min           = 0.001;
    double frequency_max           = 1.0;
    int    nfrequencies            = 0;
    double peak_enhancement_factor = 1.0;
    bool   is_normalized           = false;
    int    seed                    = 1;
    bool   wave_stretching         = true;
    double ramp_duration           = 0.0;
    std::string eta_file_path;
};

class IrregularWaves : public WaveBase {
  public:
    explicit IrregularWaves(const IrregularWaveParams& params);
    void Initialize() override {}

    void CreateSpectrum();
    std::vector<double> GetSpectrum();
    const std::vector<double>& GetFreeSurfaceElevation() const;
    const std::vector<double>& GetFreeSurfaceTime() const;
    std::vector<double> GetFrequenciesHz() const;

    std::pair<std::vector<double>, std::vector<double>>
    ComputeElevationTimeSeries(double t_start, double t_end, double dt) const;

    Eigen::VectorXd GetForceAtTime(double t) const override;
    WaveMode GetWaveMode() const override { return mode_; }

    Eigen::VectorXd SetSpectrumFrequencies(double start, double end, int num_steps);

    void AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data, const HydroData::SimulationParameters& sim_data);

    double GetIRFTimeMax() const { return irf_time_max_; }

    double GetElevation(const Eigen::Vector3d& position, double time) const override;
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time, double elevation) const override;
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) const override;

    /// Return the surface slope (∂η/∂x, ∂η/∂y) at a given position and time.
    /// Used for computing surface normals in visualization.
    Eigen::Vector2d GetElevationGradientXY(const Eigen::Vector3d& position, double time) const;

    /// Compute elevation using only the first max_components wave components.
    /// Used for faster visualization rendering when full accuracy is not required.
    /// @param max_components Maximum number of frequency components to use (0 or negative = use all)
    double GetElevationForVisualization(const Eigen::Vector3d& position, double time, int max_components) const;

  private:
    IrregularWaveParams params_;
    bool spectrum_created_ = false;
    bool use_eta_from_file_ = false;

    static constexpr WaveMode mode_ = WaveMode::irregular;
    std::vector<HydroData::IrregularWaveInfo> wave_info_;
    std::vector<Eigen::MatrixXd> ex_irf_sampled_;
    std::vector<Eigen::VectorXd> ex_irf_time_sampled_;
    std::vector<Eigen::VectorXd> ex_irf_width_sampled_;
    Eigen::VectorXd spectrum_frequencies_;
    Eigen::VectorXd spectral_densities_;
    Eigen::VectorXd spectral_widths_;
    Eigen::VectorXd wavenumbers_;
    Eigen::VectorXd wave_phases_;

    // Stored eta data — only populated when importing from file.
    std::vector<double> free_surface_elevation_sampled_;
    std::vector<double> free_surface_time_sampled_;

    // Pre-computed arrays for fast elevation calculation (SIMD-optimized).
    // Computed once from the spectrum and reused for every GetElevation() call.
    Eigen::VectorXd amplitudes_;
    Eigen::VectorXd angular_freqs_;

    // Pre-computed excitation transfer function coefficients.
    // Converts the O(N_irf * N_freq) time-domain convolution into an
    // O(N_freq) frequency-domain evaluation per DOF per timestep:
    //   F_ex(dof, t) = sum_i A_i * [C(dof,i)*cos(theta_i) - S(dof,i)*sin(theta_i)]
    // where theta_i = phi_i - omega_i * t.
    std::vector<Eigen::MatrixXd> ex_transfer_cos_;  // [body](6, N_freq)
    std::vector<Eigen::MatrixXd> ex_transfer_sin_;  // [body](6, N_freq)

    // Pre-computed trig matrices for batch eta evaluation during ramp period.
    // irf_cos_wt_[b](j, i) = cos(omega_i * tau_j), similarly for sin.
    // Allows computing eta at all IRF points via matrix-vector product.
    std::vector<Eigen::MatrixXd> irf_cos_wt_;  // [body](N_irf, N_freq)
    std::vector<Eigen::MatrixXd> irf_sin_wt_;  // [body](N_irf, N_freq)
    double irf_time_max_ = 0.0;

    void InitializeIRFVectors();
    void PrecomputeAmplitudes();
    void PrecomputeExcitationTransfer();
    std::vector<double> ReadEtaFromFile();

    const Eigen::MatrixXd& GetExcitationIRF(int b) const;
    void CalculateWidthIRF();
    double ExcitationConvolution(int body, int dof, double time) const;
};

#endif  // HYDROC_WAVES_IRREGULAR_WAVE_H

