/*********************************************************************
 * @file wave_types.cpp
 *
 * @brief implementation file for Wavebase and classes inheriting from WaveBase.
 *********************************************************************/
#include <hydroc/wave_types.h>
#include <unsupported/Eigen/Splines>
#include <hydroc/helper.h>

Eigen::VectorXd NoWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies_ * 6;
    Eigen::VectorXd f(dof);
    for (int i = 0; i < dof; i++) {
        f[i] = 0.0;
    }
    return f;
}

RegularWave::RegularWave() {
    num_bodies_ = 1;
}

RegularWave::RegularWave(unsigned int num_b) {
    num_bodies_ = num_b;
}

void RegularWave::Initialize() {
    // set up regular waves here, call other helper functions as necessary
    int total_dofs = 6 * num_bodies_;
    excitation_force_mag_.resize(total_dofs);
    excitation_force_phase_.resize(total_dofs);
    force_.resize(total_dofs);

    double wave_omega_delta = GetOmegaDelta();
    double freq_index_des   = (regular_wave_omega_ / wave_omega_delta) - 1;
    for (int b = 0; b < num_bodies_; b++) {
        for (int rowEx = 0; rowEx < 6; rowEx++) {
            int body_offset = 6 * b;
            // why are these always 0? vvv TODO check/change this
            excitation_force_mag_[body_offset + rowEx]   = GetExcitationMagInterp(b, rowEx, 0, freq_index_des);
            excitation_force_phase_[body_offset + rowEx] = GetExcitationPhaseInterp(b, rowEx, 0, freq_index_des);
        }
    }
}

void RegularWave::AddH5Data(std::vector<HydroData::RegularWaveInfo>& reg_h5_data) {
    wave_info_ = reg_h5_data;
}

Eigen::VectorXd RegularWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies_ * 6;
    Eigen::VectorXd f(dof);
    // initialize the force here:
    for (int b = 0; b < num_bodies_; b++) {
        int body_offset = 6 * b;
        for (int rowEx = 0; rowEx < 6; rowEx++) {
            f[body_offset + rowEx] = excitation_force_mag_[body_offset + rowEx] * regular_wave_amplitude_ *
                                     cos(regular_wave_omega_ * t + excitation_force_phase_[rowEx]);
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
    double freq_interp_val    = freq_index_des - floor(freq_index_des);
    double excitationMagFloor = wave_info_[b].excitation_mag_matrix(i, j, (int)floor(freq_index_des));
    double excitationMagCeil  = wave_info_[b].excitation_mag_matrix(i, j, (int)floor(freq_index_des) + 1);
    double excitationMag      = (freq_interp_val * (excitationMagCeil - excitationMagFloor)) + excitationMagFloor;

    return excitationMag;
}

double RegularWave::GetExcitationPhaseInterp(int b, int i, int j, double freq_index_des) const {
    double freq_interp_val      = freq_index_des - floor(freq_index_des);  // look into c++ modf TODO
    double excitationPhaseFloor = wave_info_[b].excitation_phase_matrix(
        i, j, (int)floor(freq_index_des));  // TODO check if freq_index_des is >0, if so just cast instead of floor
    double excitationPhaseCeil = wave_info_[b].excitation_phase_matrix(i, j, (int)floor(freq_index_des) + 1);
    double excitationPhase = (freq_interp_val * (excitationPhaseCeil - excitationPhaseFloor)) + excitationPhaseFloor;

    return excitationPhase;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//// Irregular wave class definitions:
//IrregularWave::IrregularWave() {
//    num_bodies = 1;
//}
//
//IrregularWave::IrregularWave(unsigned int num_b) {
//    num_bodies = num_b;
//}
//
//void IrregularWave::Initialize() {
//    std::vector<Eigen::MatrixXd> ex_irf_old(num_bodies);
//    std::vector<Eigen::VectorXd> ex_irf_time_old(num_bodies);
//    for (int b = 0; b < num_bodies; b++) {
//        ex_irf_old[b]      = GetExcitationIRF(b);
//        ex_irf_time_old[b] = wave_info[b].excitation_irf_time;
//    }
//
//    // resample excitation IRF time series
//    // h5 file irf has different timestep, want to resample with interpolation (cubic spline?) once at start no
//    // interpolation in convolution integral part
//    // different one for each body it's 6x1x1000 so maybe switch to 2d reading
//    ex_irf_resampled.resize(num_bodies);
//    ex_irf_time_resampled.resize(num_bodies);
//    for (int b = 0; b < num_bodies; b++) {  // function call (time_in, vals_in, &time_out, &vals_out)
//        ex_irf_time_resampled[b] = ResampleTime(ex_irf_time_old[b], simulation_dt);
//        ex_irf_resampled[b]      = ResampleVals(ex_irf_time_old[b], ex_irf_old[b], ex_irf_time_resampled[b]);
//    }
//
//    CreateSpectrum();              // output in spectral_densities.txt
//    CreateFreeSurfaceElevation();  // eta initialized in here, and output to eta.txt
//}
//
//Eigen::MatrixXd IrregularWave::ResampleVals(const Eigen::VectorXd& t_old,
//                                            Eigen::MatrixXd& vals_old,
//                                            const Eigen::VectorXd& t_new) {
//    assert(vals_old.rows() == 6);
//
//    Eigen::MatrixXd vals_new(6, t_new.size());
//    // we need to scale t to be [0,1] for spline use
//    // TODO: change this to accomodate variable dt instead of const dt
//    Eigen::VectorXd t_old_scaled = Eigen::VectorXd::LinSpaced(t_old.size(), 0, 1);
//    Eigen::VectorXd t_new_scaled = Eigen::VectorXd::LinSpaced(t_new.size(), 0, 1);
//
//    Eigen::Spline<double, 6> spline =
//        Eigen::SplineFitting<Eigen::Spline<double, 6>>::Interpolate(vals_old, 3, t_old_scaled);
//    for (int i = 0; i < t_new.rows(); i++) {
//        vals_new.col(i) = spline(t_new_scaled[i]);
//    }
//
//    return vals_new;
//}
//
//Eigen::VectorXd IrregularWave::ResampleTime(const Eigen::VectorXd& t_old, const double dt_new) {
//    double dt_old    = t_old[1] - t_old[0];
//    int size_new     = static_cast<int>(ceil(t_old.size() * dt_old / dt_new));
//    double t_initial = t_old[0];
//    double t_final   = t_old[t_old.size() - 1];
//    // t_0 and t_final should be the same in old/new versions, use new time step, which means a different number of
//    // timesteps
//    Eigen::VectorXd t_new = Eigen::VectorXd::LinSpaced(size_new, t_initial, t_final);
//
//    // print for testing
//    int N_old = t_old.size() - 1;
//    int N_new = t_new.size() - 1;
//    //std::cout << "T_old = {t_i | i = 0 .. " << N_old << ", t_0 = " << t_old[0] << ", t_" << N_old << " = "
//    //          << t_old[N_old] << "}\n"
//    //          << "dt_old = " << dt_old << "\nT_new = {t_i | i = 0 .. " << N_new << ", t_0 = " << t_new[0] << ", t_"
//    //          << N_new << " = " << t_new[N_new] << "}\n"
//    //          << "dt_new = " << dt_new << std::endl;
//
//    return t_new;
//}
//
//void IrregularWave::AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data,
//                              HydroData::SimulationParameters& sim_data) {
//    wave_info = irreg_h5_data;
//    sim_data  = sim_data;
//}
//
//Eigen::VectorXd IrregularWave::GetForceAtTime(double t) {
//    unsigned int total_dofs = num_bodies * 6;
//    Eigen::VectorXd f(total_dofs);
//    // initialize the force here:
//    // for now set f to all zeros
//    for (int i = 0; i < total_dofs; i++) {
//        f[i] = 0.0;
//    }
//    // see ComputeForceExcitation and convolution functions
//
//    // force_excitation.resize(total_dofs, 0.0);
//
//    for (int body = 0; body < num_bodies; body++) {
//        // Loop through the DOFs
//        for (int dof = 0; dof < 6; ++dof) {
//            // Compute the convolution for the current DOF
//            double f_dof          = ExcitationConvolution(body, dof, t);
//            unsigned int b_offset = body * 6;
//            f[b_offset + dof]     = f_dof;
//        }
//    }
//    return f;
//}
//
/////*******************************************************************************
//// * IrregularWave::GetExcitationIRF()
//// * returns the std::vector of excitation_irf_matrix from h5 file
//// *******************************************************************************/
////Eigen::MatrixXd IrregularWave::GetExcitationIRF(int b) const {
////    return wave_info[b].excitation_irf_matrix;
////}
//
//double IrregularWave::ExcitationConvolution(int body, int dof, double time) {
//    double f_ex  = 0.0;
//    double width = ex_irf_time_resampled[body][1] - ex_irf_time_resampled[body][0];
//    //std::cout << "width=" << width << std::endl;
//
//    for (size_t j = 0; j < ex_irf_time_resampled[0].size(); ++j) {
//        double tau        = ex_irf_time_resampled[0][j];
//        double t_tau      = time - tau;
//        double ex_irf_val = ex_irf_resampled[body](dof, j);
//        if (0.0 < t_tau && t_tau < eta.size() * width) {
//            size_t eta_index = static_cast<size_t>(t_tau / width);
//            double eta_val   = eta[eta_index - 1];
//            f_ex += ex_irf_val * eta_val * width;  // eta is wave elevation
//        }
//    }
//
//    return f_ex;
//}
//
//Eigen::VectorXd IrregularWave::SetSpectrumFrequencies(double start, double end, int num_points) {
//    Eigen::VectorXd result(num_points);
//    double step = (end - start) / (num_points - 1);
//
//    for (int i = 0; i < num_points; ++i) {
//        result[i] = start + i * step;
//    }
//
//    spectrum_frequencies = result;
//
//    return result;
//}
//
//void IrregularWave::CreateSpectrum() {
//    // Define the frequency vector
//    spectrum_frequencies = Eigen::VectorXd::LinSpaced(1000, 0.001, 1.0);
//
//    // Calculate the Pierson-Moskowitz Spectrum
//    spectral_densities = PiersonMoskowitzSpectrumHz(spectrum_frequencies, wave_height, wave_period);
//
//    // Open a file stream for writing
//    std::ofstream outputFile("spectral_densities.txt");
//
//    // Check if the file stream is open
//    if (outputFile.is_open()) {
//        // Write the spectral densities and their corresponding frequencies to the file
//        for (size_t i = 0; i < spectral_densities.size(); ++i) {
//            outputFile << spectrum_frequencies[i] << " : " << spectral_densities[i] << std::endl;
//        }
//
//        // Close the file stream
//        outputFile.close();
//    } else {
//        std::cerr << "Unable to open file for writing." << std::endl;
//    }
//}
//
//void IrregularWave::CreateFreeSurfaceElevation() {
//    // Create a time index vector
//    // UpdateNumTimesteps();
//    int num_timesteps = static_cast<int>(simulation_duration / simulation_dt) + 1;
//
//    Eigen::VectorXd time_index = Eigen::VectorXd::LinSpaced(num_timesteps, 0, simulation_duration);
//
//    // Calculate the free surface elevation
//    eta = FreeSurfaceElevation(spectrum_frequencies, spectral_densities, time_index, sim_data.water_depth);
//
//    // Apply ramp if ramp_duration is greater than 0
//    if (ramp_duration > 0.0) {
//        // UpdateRampTimesteps();
//        int ramp_timesteps   = static_cast<int>(ramp_duration / simulation_dt) + 1;
//        Eigen::VectorXd ramp = Eigen::VectorXd::LinSpaced(ramp_timesteps, 0.0, 1.0);
//
//        for (size_t i = 0; i < ramp.size(); ++i) {
//            eta[i] *= ramp[i];
//        }
//    }
//
//    // Open a file stream for writing
//    std::ofstream eta_output("eta.txt");
//    // Check if the file stream is open
//    if (eta_output.is_open()) {
//        // Write the spectral densities and their corresponding frequencies to the file
//        for (size_t i = 0; i < eta.size(); ++i) {
//            eta_output << time_index[i] << " : " << eta[i] << std::endl;
//        }
//        // Close the file stream
//        eta_output.close();
//    } else {
//        std::cerr << "Unable to open file for writing." << std::endl;
//    }
//
//    // std::vector<std::array<double, 3>> free_surface_3d_pts    = CreateFreeSurface3DPts(eta, time_index);
//    // std::vector<std::array<size_t, 3>> free_surface_triangles = CreateFreeSurfaceTriangles(time_index.size());
//
//    // WriteFreeSurfaceMeshObj(free_surface_3d_pts, free_surface_triangles, "fse_mesh.obj");
//}


std::vector<double> ComputeWaveNumbers(const std::vector<double>& omegas,
                                       double water_depth,
                                       double g           = 9.81,
                                       double tolerance   = 1e-6,
                                       int max_iterations = 100) {
    std::vector<double> wave_numbers(omegas.size());

    for (size_t i = 0; i < omegas.size(); ++i) {
        double omega = omegas[i];

        // Initial guess for wave number (using deep water approximation)
        double k = omega * omega / g;

        int iterations = 0;
        double error   = 1.0;
        while (error > tolerance && iterations < max_iterations) {
            double tanh_kh = std::tanh(k * water_depth);
            double f       = omega * omega - g * k * tanh_kh;
            double df      = -2.0 * g * tanh_kh - g * k * water_depth * (1.0 - tanh_kh * tanh_kh);

            double delta_k = f / df;
            k -= delta_k;
            error = std::abs(delta_k);
            iterations++;
        }

        wave_numbers[i] = k;
    }

    return wave_numbers;
}

std::vector<double> FreeSurfaceElevation(const Eigen::VectorXd& freqs_hz,
                                     const Eigen::VectorXd& spectral_densities,
                                     const Eigen::VectorXd& time_index,
                                     double water_depth,
                                     int seed) {
    double delta_f = freqs_hz(Eigen::last) / freqs_hz.size();
    std::vector<double> omegas(freqs_hz.size());

    for (size_t i = 0; i < freqs_hz.size(); ++i) {
        omegas[i] = 2 * M_PI * freqs_hz[i];
    }

    std::vector<double> wave_numbers = ComputeWaveNumbers(omegas, water_depth);

    std::vector<double> A(spectral_densities.size());
    for (size_t i = 0; i < spectral_densities.size(); ++i) {
        A[i] = 2 * spectral_densities[i] * delta_f;
    }

    std::vector<double> sqrt_A(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        sqrt_A[i] = std::sqrt(A[i]);
    }

    std::vector<std::vector<double>> omegas_t(time_index.size(), std::vector<double>(omegas.size()));
    for (size_t i = 0; i < time_index.size(); ++i) {
        for (size_t j = 0; j < omegas.size(); ++j) {
            omegas_t[i][j] = time_index[i] * omegas[j];
        }
    }

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 2 * M_PI);
    std::vector<double> phases(omegas.size());
    for (size_t i = 0; i < phases.size(); ++i) {
        phases[i] = dist(rng);
    }

    std::vector<double> eta(time_index.size(), 0.0);
    for (size_t i = 0; i < spectral_densities.size(); ++i) {
        for (size_t j = 0; j < time_index.size(); ++j) {
            eta[j] += sqrt_A[i] * std::cos(omegas_t[j][i] + phases[i]);
        }
    }

    return eta;
}

//void IrregularWave::SetUpWaveMesh(std::string filename) {
//    mesh_file_name             = filename;
//    int num_timesteps          = static_cast<int>(simulation_duration / simulation_dt) + 1;
//    Eigen::VectorXd time_index = Eigen::VectorXd::LinSpaced(num_timesteps, 0, simulation_duration);
//    std::vector<std::array<double, 3>> free_surface_3d_pts    = CreateFreeSurface3DPts(eta, time_index);
//    std::vector<std::array<size_t, 3>> free_surface_triangles = CreateFreeSurfaceTriangles(time_index.size());
//
//    WriteFreeSurfaceMeshObj(free_surface_3d_pts, free_surface_triangles, mesh_file_name);
//}

std::vector<std::array<double, 3>> CreateFreeSurface3DPts(const std::vector<double>& eta,
                                                          const Eigen::VectorXd& t_vec) {
    std::vector<std::array<double, 3>> surface(t_vec.size() * 2);

    for (size_t i = 0; i < t_vec.size(); ++i) {
        double t = -1 * t_vec[i];
        double z = eta[i];

        surface[2 * i]     = {t, -10.0, z};
        surface[2 * i + 1] = {t, 10.0, z};
    }

    return surface;
}

std::vector<std::array<size_t, 3>> CreateFreeSurfaceTriangles(size_t eta_size) {
    std::vector<std::array<size_t, 3>> triangles;

    for (size_t i = 0; i < eta_size / 2 - 1; ++i) {
        triangles.push_back({2 * i, 2 * i + 1, 2 * i + 3});
        triangles.push_back({2 * i, 2 * i + 3, 2 * i + 2});
    }

    return triangles;
}
#include <iomanip>

void WriteFreeSurfaceMeshObj(const std::vector<std::array<double, 3>>& points,
                             const std::vector<std::array<size_t, 3>>& triangles,
                             const std::string& file_name) {
    std::ofstream out(file_name);
    if (!out) {
        std::cerr << "Failed to open " << file_name << std::endl;
        return;
    }

    // Write header
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    out << "# Wavefront OBJ file exported by HydroChrono" << std::endl;
    out << "# File Created: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl << std::endl;

    // Write vertices
    out << "# Vertices: " << points.size() << std::endl << std::endl;
    out << std::fixed << std::setprecision(6);
    for (const auto& point : points) {
        out << "v ";
        out << std::setw(14) << point[0] << ' ';
        out << std::setw(14) << point[1] << ' ';
        out << std::setw(14) << point[2] << std::endl;
    }
    out << std::endl;

    // Write faces
    out << "# Faces: " << triangles.size() << std::endl << std::endl;
    for (const auto& triangle : triangles) {
        out << "f ";
        out << std::setw(9) << triangle[0] + 1;
        out << std::setw(9) << triangle[1] + 1;
        out << std::setw(9) << triangle[2] + 1 << std::endl;
    }

    out.close();
}

//std::string IrregularWave::GetMeshFile() {
//    return mesh_file_name;
//}
//
//Eigen::Vector3<double> IrregularWave::GetWaveMeshVelocity() {
//    return Eigen::Vector3d(1.0, 0, 0);
//}


//IrregularWaves::IrregularWaves(const std::string& eta_file_path, double wave_height, double wave_period)
//    : eta_file_path(eta_file_path), wave_height(wave_height), wave_period(wave_period) {
//    if (!eta_file_path.empty()) {
//        ReadEtaFromFile();
//        spectrumCreated = false;
//    } else {
//        CreateSpectrum();
//        CreateFreeSurfaceElevation();
//        spectrumCreated = true;
//    }
//}

IrregularWaves::IrregularWaves(const IrregularWaveParams& params)
    : num_bodies_(params.num_bodies_),
      eta_file_path_(params.eta_file_path_),
      wave_height_(params.wave_height_),
      wave_period_(params.wave_period_),
      peak_enhancement_factor_(params.peak_enhancement_factor_),
      simulation_dt_(params.simulation_dt_),
      simulation_duration_(params.simulation_duration_),
      ramp_duration_(params.ramp_duration_) {
    std::cout << "Creating IrregularWaves object..." << std::endl;
    ex_irf_resampled_.resize(num_bodies_);
    ex_irf_time_resampled_.resize(num_bodies_);
    std::cout << "eta_file_path (1)" << eta_file_path_ << std::endl;
    if (!eta_file_path_.empty()) {
        std::cout << "Reading eta.txt..." << std::endl;
        ReadEtaFromFile();
        spectrumCreated_ = false;
    } else if (wave_height_ != 0.0 && wave_period_ != 0.0) {
        CreateSpectrum();
        CreateFreeSurfaceElevation();
        spectrumCreated_ = true;
    }
    std::cout << "eta.txt read." << std::endl;

    std::vector<Eigen::MatrixXd> ex_irf_old(num_bodies_);
    std::vector<Eigen::VectorXd> ex_irf_time_old(num_bodies_);

    //for (unsigned int b = 0; b < num_bodies; b++) {
    //    std::cout << "get excitation IRFs..." << std::endl;
    //    ex_irf_old[b]      = GetExcitationIRF(b);
    //    std::cout << "get excitation IRF time..." << std::endl;
    //    ex_irf_time_old[b] = wave_info[b].excitation_irf_time;
    //}

    //// Resample excitation IRF time series
    //std::cout << "resample IRFs..." << std::endl;
    //for (unsigned int b = 0; b < num_bodies; b++) {
    //    ex_irf_time_resampled[b] = ResampleTime(ex_irf_time_old[b], simulation_dt);
    //    ex_irf_resampled[b]      = ResampleVals(ex_irf_time_old[b], ex_irf_old[b], ex_irf_time_resampled[b]);
    //}

    std::cout << "IrregularWaves instantiated..." << std::endl;
}

void IrregularWaves::InitializeIRFVectors() {
    std::vector<Eigen::MatrixXd> ex_irf_old(num_bodies_);
    std::vector<Eigen::VectorXd> ex_irf_time_old(num_bodies_);

    for (unsigned int b = 0; b < num_bodies_; b++) {
        ex_irf_old[b]      = GetExcitationIRF(b);
        ex_irf_time_old[b] = wave_info_[b].excitation_irf_time;
    }

    // Resample excitation IRF time series
    for (unsigned int b = 0; b < num_bodies_; b++) {
        ex_irf_time_resampled_[b] = ResampleTime(ex_irf_time_old[b], simulation_dt_);
        ex_irf_resampled_[b]      = ResampleVals(ex_irf_time_old[b], ex_irf_old[b], ex_irf_time_resampled_[b]);
    }
}

std::vector<double> IrregularWaves::GetSpectrum() {
    if (!spectrumCreated_) {
        throw std::runtime_error(
            "Spectrum has not been created. Initialize with wave height and period to create spectrum.");
    }
    return spectrum_;
}

std::vector<double> IrregularWaves::GetFreeSurfaceElevation() {
    return free_surface_elevation_;
}

std::vector<double> IrregularWaves::GetEtaTimeData() {
    return time_data_;
}

void IrregularWaves::ReadEtaFromFile() {
    std::cout << "ReadEtaFromFile()" << std::endl;
    std::ifstream file(eta_file_path_);
    if (!file) {
        throw std::runtime_error("Unable to open file at: " + eta_file_path_);
    }

    std::string line;
    double time, eta;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        char delimiter;
        if (!(ss >> time >> delimiter >> eta) || delimiter != ':') {
            throw std::runtime_error("Could not parse line: " + line);
        }
        time_data_.push_back(time);
        free_surface_elevation_.push_back(eta);
    }
}

/*******************************************************************************
 * IrregularWave::GetExcitationIRF()
 * returns the std::vector of excitation_irf_matrix from h5 file
 *******************************************************************************/
Eigen::MatrixXd IrregularWaves::GetExcitationIRF(int b) const {
    return wave_info_[b].excitation_irf_matrix;
}

void IrregularWaves::AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data,
                               HydroData::SimulationParameters& sim_data) {
    wave_info_ = irreg_h5_data;
    sim_data  = sim_data;  // Note: This line seems to be redundant as mentioned before.

    //// Add this loop to print the values of wave_info[b].excitation_irf_matrix for each body
    //for (unsigned int b = 0; b < num_bodies; ++b) {
    //    std::cout << "wave_info[" << b << "].excitation_irf_matrix:\n";
    //    std::cout << wave_info[b].excitation_irf_matrix << std::endl;
    //}

    InitializeIRFVectors();
}

Eigen::VectorXd IrregularWaves::GetForceAtTime(double t) {
    unsigned int total_dofs = num_bodies_ * 6;
    Eigen::VectorXd f(total_dofs);
    // initialize the force here:
    // for now set f to all zeros
    for (int i = 0; i < total_dofs; i++) {
        f[i] = 0.0;
    }
    // see ComputeForceExcitation and convolution functions

    // force_excitation.resize(total_dofs, 0.0);

    for (int body = 0; body < num_bodies_; body++) {
        // Loop through the DOFs
        for (int dof = 0; dof < 6; ++dof) {
            // Compute the convolution for the current DOF
            double f_dof          = ExcitationConvolution(body, dof, t);
            unsigned int b_offset = body * 6;
            f[b_offset + dof]     = f_dof;
        }
    }
    
    return f;
}

Eigen::MatrixXd IrregularWaves::ResampleVals(const Eigen::VectorXd& t_old,
                                            Eigen::MatrixXd& vals_old,
                                            const Eigen::VectorXd& t_new) {
    assert(vals_old.rows() == 6);

    Eigen::MatrixXd vals_new(6, t_new.size());
    // we need to scale t to be [0,1] for spline use
    // TODO: change this to accomodate variable dt instead of const dt
    Eigen::VectorXd t_old_scaled = Eigen::VectorXd::LinSpaced(t_old.size(), 0, 1);
    Eigen::VectorXd t_new_scaled = Eigen::VectorXd::LinSpaced(t_new.size(), 0, 1);

    Eigen::Spline<double, 6> spline =
        Eigen::SplineFitting<Eigen::Spline<double, 6>>::Interpolate(vals_old, 3, t_old_scaled);
    for (int i = 0; i < t_new.rows(); i++) {
        vals_new.col(i) = spline(t_new_scaled[i]);
    }

    //std::cout << "Old time vector:\n" << t_old << std::endl;
    //std::cout << "Old values:\n" << vals_old << std::endl;
    //std::cout << "New time vector:\n" << t_new << std::endl;
    //std::cout << "New values:\n" << vals_new << std::endl;  // assuming new_vals is the output of the function


    return vals_new;
}

Eigen::VectorXd IrregularWaves::ResampleTime(const Eigen::VectorXd& t_old, const double dt_new) {
    double dt_old    = t_old[1] - t_old[0];
    int size_new     = static_cast<int>(ceil(t_old.size() * dt_old / dt_new));
    double t_initial = t_old[0];
    double t_final   = t_old[t_old.size() - 1];
    // t_0 and t_final should be the same in old/new versions, use new time step, which means a different number of
    // timesteps
    Eigen::VectorXd t_new = Eigen::VectorXd::LinSpaced(size_new, t_initial, t_final);

    int N_old = t_old.size() - 1;
    int N_new = t_new.size() - 1;

    return t_new;
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
    // Define the frequency vector
    spectrum_frequencies_ = Eigen::VectorXd::LinSpaced(1000, 0.001, 1.0);

    // Calculate the Pierson-Moskowitz Spectrum
    spectral_densities_ = JONSWAPSpectrumHz(spectrum_frequencies_, wave_height_, wave_period_, peak_enhancement_factor_);

    // Open a file stream for writing
    std::ofstream outputFile("spectral_densities.txt");

    // Check if the file stream is open
    if (outputFile.is_open()) {
        // Write the spectral densities and their corresponding frequencies to the file
        for (size_t i = 0; i < spectral_densities_.size(); ++i) {
            outputFile << spectrum_frequencies_[i] << " : " << spectral_densities_[i] << std::endl;
        }

        // Close the file stream
        outputFile.close();
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}

// TODO put spectrum functions in a new namespace (when we have more options?)
Eigen::VectorXd PiersonMoskowitzSpectrumHz(Eigen::VectorXd& f, double Hs, double Tp) {
    // Sort the frequency vector
    std::sort(f.begin(), f.end());

    // Initialize the spectral densities vector
    Eigen::VectorXd spectral_densities(f.size());

    // Calculate the spectral densities
    for (size_t i = 0; i < f.size(); ++i) {
        spectral_densities[i] = 1.25 * std::pow(1 / Tp, 4) * std::pow(Hs / 2, 2) * std::pow(f[i], -5) *
                                std::exp(-1.25 * std::pow(1 / Tp, 4) * std::pow(f[i], -4));
    }

    return spectral_densities;
}

Eigen::VectorXd JONSWAPSpectrumHz(Eigen::VectorXd& f, double Hs, double Tp, double gamma) {
    auto spectral_densities = PiersonMoskowitzSpectrumHz(f, Hs, Tp);

    // Scale spectral densities from PM to JONSWAP with gamma factor
    for (size_t i = 0; i < spectral_densities.size(); ++i) {
        double sigma;
        if (f[i] <= 1.0 / Tp) {
            sigma = 0.07;
        } else {
            sigma = 0.09;
        }
        spectral_densities[i] *= pow(gamma, exp(-(1.0 / (2.0 * pow(sigma, 2))) * pow(f[i] * Tp - 1.0, 2)));
    }

    return spectral_densities;
}

void IrregularWaves::CreateFreeSurfaceElevation() {
    // Create a time index vector
    // UpdateNumTimesteps();
    int num_timesteps = static_cast<int>(simulation_duration_ / simulation_dt_) + 1;

    Eigen::VectorXd time_index = Eigen::VectorXd::LinSpaced(num_timesteps, 0, simulation_duration_);

    // Calculate the free surface elevation
    free_surface_elevation_ = FreeSurfaceElevation(spectrum_frequencies_, spectral_densities_, time_index, sim_data_.water_depth);

    // Apply ramp if ramp_duration is greater than 0
    if (ramp_duration_ > 0.0) {
        // UpdateRampTimesteps();
        int ramp_timesteps   = static_cast<int>(ramp_duration_ / simulation_dt_) + 1;
        Eigen::VectorXd ramp = Eigen::VectorXd::LinSpaced(ramp_timesteps, 0.0, 1.0);

        for (size_t i = 0; i < ramp.size(); ++i) {
            free_surface_elevation_[i] *= ramp[i];
        }
    }

    // Open a file stream for writing
    std::ofstream eta_output("eta.txt");
    // Check if the file stream is open
    if (eta_output.is_open()) {
        // Write the spectral densities and their corresponding frequencies to the file
        for (size_t i = 0; i < free_surface_elevation_.size(); ++i) {
            eta_output << time_index[i] << " : " << free_surface_elevation_[i] << std::endl;
        }
        // Close the file stream
        eta_output.close();
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }

    // std::vector<std::array<double, 3>> free_surface_3d_pts    = CreateFreeSurface3DPts(eta, time_index);
    // std::vector<std::array<size_t, 3>> free_surface_triangles = CreateFreeSurfaceTriangles(time_index.size());

    // WriteFreeSurfaceMeshObj(free_surface_3d_pts, free_surface_triangles, "fse_mesh.obj");
}

double IrregularWaves::ExcitationConvolution(int body, int dof, double time) {
    double f_ex  = 0.0;
    double width = ex_irf_time_resampled_[body][1] - ex_irf_time_resampled_[body][0];
    // std::cout << "width=" << width << std::endl;
    for (size_t j = 0; j < ex_irf_time_resampled_[0].size(); ++j) {
        double tau        = ex_irf_time_resampled_[0][j];
        double t_tau      = time - tau;
        
        double ex_irf_val = ex_irf_resampled_[body](dof, j);

        if (0.0 < t_tau && t_tau < free_surface_elevation_.size() * width) {
            size_t eta_index = static_cast<size_t>(t_tau / width);
            double eta_val   = free_surface_elevation_[eta_index - 1];

            //std::cout << "ex_irf_val = " << ex_irf_val << std::endl;
            //std::cout << "eta_val = " << eta_val << std::endl;
            //std::cout << "width = " << width << std::endl;

            f_ex += ex_irf_val * eta_val * width;  // eta is wave elevation
        }
    }

    return f_ex;
}

 void IrregularWaves::SetUpWaveMesh(std::string filename) {
    mesh_file_name_             = filename;
    int num_timesteps          = static_cast<int>(simulation_duration_ / simulation_dt_) + 1;
    Eigen::VectorXd time_index = Eigen::VectorXd::LinSpaced(num_timesteps, 0, simulation_duration_);
    std::vector<std::array<double, 3>> free_surface_3d_pts = CreateFreeSurface3DPts(free_surface_elevation_, time_index);
    std::vector<std::array<size_t, 3>> free_surface_triangles = CreateFreeSurfaceTriangles(time_index.size());

    WriteFreeSurfaceMeshObj(free_surface_3d_pts, free_surface_triangles, mesh_file_name_);
}

std::string IrregularWaves::GetMeshFile() {
    return mesh_file_name_;
}

Eigen::Vector3<double> IrregularWaves::GetWaveMeshVelocity() {
    return Eigen::Vector3d(1.0, 0, 0);
}