#pragma once

#include <optional>
#include <string>
#include <vector>
#include <variant>
#include <limits>

namespace H5 {
class H5File;
}

#include <chrono/core/ChMatrix.h>
#include <unsupported/Eigen/CXX11/Tensor>

/**@brief Extract bemio formated hdf5 data
 *
 * https://wec-sim.github.io/bemio/_modules/bemio/io/output.html
 */

class H5FileInfo;

// contains "chunked" data from the h5 file, pass this class to the h5fileinfo class to populate entries
class HydroData {
  public:
    struct BodyInfo {
        std::string body_name;
        int body_num;
        double disp_vol;
        Eigen::VectorXd rirf_time_vector;
        double rirf_timestep;
        Eigen::VectorXd cg;
        Eigen::VectorXd cb;
        Eigen::MatrixXd lin_matrix;
        Eigen::MatrixXd inf_added_mass;
        Eigen::Tensor<double, 3> rirf_matrix;
        // Eigen::Tensor<double, 3> radiation_damping_matrix;
    };
    struct SimulationParameters {
        std::string h5_file_name;
        double rho;
        double g;
        double water_depth;
    };
    struct RegularWaveInfo {
        Eigen::VectorXd freq_list;
        Eigen::Tensor<double, 3> excitation_mag_matrix;
        Eigen::Tensor<double, 3> excitation_phase_matrix;
    };
    struct IrregularWaveInfo {
        // Eigen::Tensor<double,3> excitation_re_matrix;
        // Eigen::Vector3i re_dims;
        // Eigen::Tensor<double, 3> excitation_im_matrix;
        // Eigen::Vector3i im_dims;
        Eigen::VectorXd excitation_irf_time;
        Eigen::MatrixXd excitation_irf_matrix;  // TODO needs to be tensor?

        // see std::optional documentation for how to use
        std::optional<Eigen::MatrixXd> excitation_irf_resampled;  // TODO needs to be tensor?
        std::optional<Eigen::MatrixXd> excitation_irf_time_resampled;
    };

  private:
    // a vector of BodyInfo objects, one for each hydro body in system
    std::vector<BodyInfo> body_data;
    // a variable for simulation parameters, assumed the same for each body, only  1 per system
    SimulationParameters sim_data;
    // a vector of RegularWaveInfo, one for each hydro body in system
    // is empty if regular waves are not used
    std::vector<RegularWaveInfo> reg_wave_data;
    // a vector of IrregularWaveInfo, one for each hydro body in system
    // is empty if irregular waves are not used
    std::vector<IrregularWaveInfo> irreg_wave_data;
    friend H5FileInfo;
    void resize(int num_bodies);
    HydroData() = default;

  public:
    // getter function naming conventions: put Matrix, Vector, or Val at end to denote return type, always put body_num
    // argument first, body_num is always 0 indexed, one line return types can be defined here and not in cpp
    Eigen::MatrixXd GetInfAddedMassMatrix(int b) const;
    double GetHydrostaticStiffnessVal(int b, int i, int j) const;
    Eigen::MatrixXd GetLinMatrix(int b) const;
    double GetRIRFVal(int b, int i, int n, int m) const;
    double GetDispVolVal(int b) const { return body_data[b].disp_vol; }
    Eigen::VectorXd GetCGVector(int b) const { return body_data[b].cg; }
    Eigen::VectorXd GetCBVector(int b) const { return body_data[b].cb; }
    double GetExcitationIRFVal(int b, int dof, int s) const;
    Eigen::MatrixXd GetExcitationIRF(int b) const;
    //Eigen::VectorXd ResampleExcitationIRFTime(double dt_new);
    //Eigen::MatrixXd GetExcitationIRFResampled(int b) const;
    //std::pair<Eigen::VectorXd, Eigen::VectorXd> ResampleExcitationIRF(int b, double dt_new);

    // things that are the same no matter the body, don't need body argument
    int GetRIRFDims(int i) const;
    Eigen::VectorXd GetRIRFTimeVector() const;  // TODO
    double GetRhoVal() const { return sim_data.rho; }
    //Eigen::VectorXd HydroData::GetExcitationIRFTime() const;

    // getters for individual structs
    std::vector<BodyInfo>& GetBodyInfos() { return body_data; }
    SimulationParameters& GetSimulationInfo() { return sim_data; }
    std::vector<RegularWaveInfo>& GetRegularWaveInfos() { return reg_wave_data; }
    std::vector<IrregularWaveInfo>& GetIrregularWaveInfos() { return irreg_wave_data; }
};

class H5FileInfo {
  public:
    bool printed = false;

    H5FileInfo(std::string file, int num_bod = 1);
    H5FileInfo() = delete;

    H5FileInfo(const H5FileInfo& old) = default;
    H5FileInfo& operator=(const H5FileInfo& rhs) = default;

    H5FileInfo(H5FileInfo&&) = default;
    H5FileInfo& operator=(H5FileInfo&& rhs) = default;

    ~H5FileInfo();

    HydroData readH5Data();  // TODO: eventually pass user input struct here

  private:
    std::string h5_file_name;
    int num_bodies;

    //void InitScalar(H5::H5File& file, std::string data_name, double& var);
    void InitScalar(H5::H5File& file, std::string data_name, double& var);
    void Init1D(H5::H5File& file, std::string data_name, Eigen::VectorXd& var);
    void Init2D(H5::H5File& file, std::string data_name, Eigen::MatrixXd& var);
    void Init3D(H5::H5File& file, std::string data_name, Eigen::Tensor<double, 3>& var /*, std::vector<int>& dims*/);
    Eigen::MatrixXd squeeze_mid(Eigen::Tensor<double, 3> to_be_squeezed);
};