#include <H5Cpp.h>
#include <hydroc/h5fileinfo.h>

#include <filesystem>

#include <unsupported/Eigen/Splines>

using namespace chrono;

// =============================================================================
// Misc
// =============================================================================

template <typename T>
void WriteDataToFile(const std::vector<T>& data, const std::string& filename) {
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        for (const auto& item : data) {
            outFile << item << std::endl;
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open the file for writing: " << filename << std::endl;
    }
}


// =============================================================================
// H5FileInfo Class Definitions
// =============================================================================

/*******************************************************************************
 * H5FileInfo constructor
 * requires file name (in absolute file name or referenced from executable location)
 * and body name (name of body's section in H5 file, ie "body1" etc)
 * each body in system should have its own H5FileInfo object
 * calls readH5Data()
 *******************************************************************************/
H5FileInfo::H5FileInfo(std::string file, int num_bod) {
    h5_file_name = file;
    num_bodies   = num_bod;
    std::cout << "searching for file: " << file << std::endl;
    if (std::filesystem::exists(file)) {
        std::cout << "found file at: " << std::filesystem::absolute(file) << std::endl;
    } else {
        std::cout << "h5 file does not exist, absolute file location: " << std::filesystem::absolute(file) << std::endl;
    }
}

/*******************************************************************************
 * H5FileInfo::readH5Data()
 * private member function called from constructor
 * calls Initialize functions to read h5 file information into  member variables
 *******************************************************************************/
HydroData H5FileInfo::readH5Data() {
    // open file with read only access
    H5::H5File userH5File(h5_file_name, H5F_ACC_RDONLY);
    HydroData data_to_init;
    data_to_init.resize(num_bodies);

    // simparams first
    InitScalar(userH5File, "simulation_parameters/rho", data_to_init.sim_data.rho);
    InitScalar(userH5File, "simulation_parameters/g", data_to_init.sim_data.g);

    // for each body things
    for (int i = 0; i < num_bodies; i++) {
        // body data
        data_to_init.body_data[i].body_name = "body" + std::to_string(i + 1);
        std::string bodyName                = data_to_init.body_data[i].body_name;  // shortcut for reading later
        data_to_init.body_data[i].body_num  = i;

        InitScalar(userH5File, bodyName + "/properties/disp_vol", data_to_init.body_data[i].disp_vol);
        Init1D(userH5File, bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/t",
               data_to_init.body_data[i].rirf_time_vector);

        // TODO add these in after the merge
    // Init1D(userH5File, bodyName + "/hydro_coeffs/excitation/impulse_response_fun/t", excitation_irf_time);
    // Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/impulse_response_fun/f", excitation_irf_matrix,
    //       excitation_irf_dims);

        // do not need rirf_timestep?
        data_to_init.body_data[i].rirf_timestep =
            data_to_init.body_data[i].rirf_time_vector[1] - data_to_init.body_data[i].rirf_time_vector[0];

        Init1D(userH5File, bodyName + "/properties/cg", data_to_init.body_data[i].cg);
        Init1D(userH5File, bodyName + "/properties/cb", data_to_init.body_data[i].cb);
        Init2D(userH5File, bodyName + "/hydro_coeffs/linear_restoring_stiffness", data_to_init.body_data[i].lin_matrix);
        Init2D(userH5File, bodyName + "/hydro_coeffs/added_mass/inf_freq", data_to_init.body_data[i].inf_added_mass);
        Init3D(userH5File, bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/K",
               data_to_init.body_data[i].rirf_matrix);
        // Init3D(userH5File, bodyName + "/hydro_coeffs/radiation_damping/all",
        //       data_to_init.body_data[i].radiation_damping_matrix);

        // reg wave
        Init1D(userH5File, "simulation_parameters/w", data_to_init.reg_wave_data[i].freq_list);
        Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/mag",
               data_to_init.reg_wave_data[i].excitation_mag_matrix);
        Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/phase",
               data_to_init.reg_wave_data[i].excitation_phase_matrix);

        // irreg wave
        // Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/re", excitation_re_matrix, re_dims);
        // Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/im", excitation_im_matrix, im_dims);
    }

    userH5File.close();
    // WriteDataToFile(excitation_irf_dims, "excitation_irf_dims.txt");
    // WriteDataToFile(excitation_irf_matrix, "excitation_irf_matrix.txt");
    return data_to_init;
}

/*******************************************************************************
 * H5FileInfo::InitScalar(file, data_name, var)
 * file - H5File& to open h5 file with body information
 * data_name - string name of H5::DataSet to be read within h5 file
 * var - member variable to store the information in data_name in
 * Reads a double type variable from data_name DataSet in file, stores it in var
 *******************************************************************************/
void H5FileInfo::InitScalar(H5::H5File& file, std::string data_name, double& var) {
    H5::DataSet dataset     = file.openDataSet(data_name);
    H5::DataSpace filespace = dataset.getSpace();
    hsize_t dims[2]         = {0, 0};
    int rank                = filespace.getSimpleExtentDims(dims);
    H5::DataSpace mspace1   = H5::DataSpace(rank, dims);
    double* temp            = new double[dims[0] * dims[1]];
    dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
    var = temp[0];
    dataset.close();
    delete[] temp;
}

/*******************************************************************************
 * H5FileInfo::Init1D(file, data_name, var)
 * file - H5File& to open h5 file with body information
 * data_name - string name of H5::DataSet to be read within h5 file
 * var - member variable to store the information in data_name in
 * Reads a 1D double type variable from data_name DataSet in file, stores it in var
 *******************************************************************************/
void H5FileInfo::Init1D(H5::H5File& file, std::string data_name, Eigen::VectorXd& var) {
    // open specific dataset
    H5::DataSet dataset = file.openDataSet(data_name);
    // Get filespace for rank and dimension
    H5::DataSpace filespace = dataset.getSpace();
    // Get number of dimensions in the file dataspace
    // Get and print the dimension sizes of the file dataspace
    hsize_t dims[2] = {0, 0};  // dataset dimensions
    int rank        = filespace.getSimpleExtentDims(dims);
    // read file into data_out 2d array
    H5::DataSpace mspace1(rank, dims);
    double* temp = new double[dims[0] * dims[1]];
    // read file info into current_pos
    dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
    var.resize(dims[0] * dims[1]);
    for (int i = 0; i < dims[0] * dims[1]; i++) {
        var[i] = temp[i];
    }
    dataset.close();
    delete[] temp;
}

/*******************************************************************************
 * H5FileInfo::Init2D(file, data_name, var)
 * file - H5File& to open h5 file with body information
 * data_name - string name of H5::DataSet to be read within h5 file
 * var - member variable to store the information in data_name in
 * Reads a 2D double type variable from data_name DataSet in file, stores it in var
 *******************************************************************************/
void H5FileInfo::Init2D(H5::H5File& file, std::string data_name, Eigen::MatrixXd& var) {
    // data_name = bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/K";
    H5::DataSet dataset     = file.openDataSet(data_name);
    H5::DataSpace filespace = dataset.getSpace();
    hsize_t dims[2]         = {0, 0};
    int rank                = filespace.getSimpleExtentDims(dims);
    // read file into data_out 2d array
    H5::DataSpace mspace(rank, dims);
    // rirf_dims[0] is number of rows, rirf_dims[1] is number of columns, rirf_dims[2] is number of matrices
    double* temp = new double[dims[0] * dims[1]];
    // read file info into data_out, a 2d array
    dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace, filespace);
    // set var here
    var.resize(dims[0], dims[1]);
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            var(i, j) = temp[i * dims[1] + j];
        }
    }
    dataset.close();
    delete[] temp;
}

/*******************************************************************************
 * H5FileInfo::Init3D(file, data_name, var)
 * file - H5File& to open h5 file with body information
 * data_name - string name of H5::DataSet to be read within h5 file
 * var - member variable to store the information in data_name in
 * d - member variable to store dimensions of 3d vectorized matrix object
 * Reads a 3D double type variable from data_name DataSet in file, stores it in var
 *******************************************************************************/
void H5FileInfo::Init3D(H5::H5File& file, std::string data_name, Eigen::Tensor<double, 3>& var) {
    // open specific dataset
    H5::DataSet dataset = file.openDataSet(data_name);
    // Get filespace for rank and dimension
    H5::DataSpace filespace = dataset.getSpace();
    hsize_t dims[3]         = {0, 0, 0};
    int rank                = filespace.getSimpleExtentDims(dims);
    // read file into data_out 2d array
    H5::DataSpace mspace(rank, dims);
    // rirf_dims[0] is number of rows, rirf_dims[1] is number of columns, rirf_dims[2] is number of matrices
    double* temp = new double[dims[0] * dims[1] * dims[2]];
    // read file info into data_out, a 2d array
    dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace, filespace);
    // set var here
    var.resize((int64_t)dims[0], (int64_t)dims[1], (int64_t)dims[2]);
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            for (int k = 0; k < dims[2]; k++) {
                int index    = k + dims[2] * (j + i * dims[1]);
                var(i, j, k) = temp[index];
            }
        }
    }
    dataset.close();
    delete[] temp;
}

/*******************************************************************************
 * H5FileInfo::~H5FileInfo()
 * H5FileInfo destructor TODO
 *******************************************************************************/
H5FileInfo::~H5FileInfo() {}

// =============================================================================
// HydroData Class Definitions
// =============================================================================

void HydroData::resize(int num_bodies) {
    body_data.resize(num_bodies);
    reg_wave_data.resize(num_bodies);
    irreg_wave_data.resize(num_bodies);
}

/*******************************************************************************
 * H5FileInfo::GetExcitationIRFDims(int i) returns the i-th component of the dimensions of excitation_irf_matrix
 * i = [0,1,2] -> [number of rows, number of columns, number of matrices]
 *******************************************************************************/
// TODO update this in next merge
// int H5FileInfo::GetExcitationIRFDims(int i) const {
//     return excitation_irf_dims[i];
// }

/*******************************************************************************
 * HydroData::GetInfAddedMassMatrix()
 * returns the matrix for added mass at infinite frequency scaled by rho for the given body
 *******************************************************************************/
Eigen::MatrixXd HydroData::GetInfAddedMassMatrix(int b) const {
    return body_data[b].inf_added_mass * sim_data.rho;
}

/*******************************************************************************
 * HydroData::GetHydrostaticStiffness()
 * returns the linear restoring stiffness matrix element in body b, row i , column j
 *******************************************************************************/
double HydroData::GetHydrostaticStiffnessVal(int b, int i, int j) const {
    return body_data[b].lin_matrix(i, j) * sim_data.rho * sim_data.g;
}

/*******************************************************************************
 * HydroData::GetLinMatrix()
 * returns the linear restoring stiffness matrix for body b
 *******************************************************************************/
Eigen::MatrixXd HydroData::GetLinMatrix(int b) const {
    return body_data[b].lin_matrix;
}

/*******************************************************************************
 * HydroData::GetRIRFval()
 * returns rirf val for DoF: 0,...,5; col: 0,...,6N-1; s: 0,...,1001 rirfdims[2]
 *******************************************************************************/
double HydroData::GetRIRFVal(int b, int dof, int col, int s) const {
    return body_data[b].rirf_matrix(dof, col, s) * sim_data.rho;  // scale radiation force by rho
}

/*******************************************************************************
 * HydroData::GetRIRFDims( int i) returns the i-th component of the dimensions of radiation_damping_matrix
 * i = [0,1,2] -> [number of rows, number of columns, number of matrices]
 *******************************************************************************/
int HydroData::GetRIRFDims(int i) const {
    return body_data[0].rirf_matrix.dimension(i);
}

/*******************************************************************************
 * HydroData::GetRIRFTimeVector()
 * returns the std::vector of rirf_time_vector from h5 file
 *******************************************************************************/
Eigen::VectorXd HydroData::GetRIRFTimeVector() const {
    return body_data[0].rirf_time_vector;
}

/*******************************************************************************
 * HydroData::GetOmegaDelta()
 * returns omega step size
 *******************************************************************************/
double HydroData::GetOmegaDelta() const {
    return GetOmegaMax() / GetNumFreqs();
}

/*******************************************************************************
 * HydroData::GetOmegaMax()
 * returns max value of omega
 *******************************************************************************/
double HydroData::GetOmegaMax() const {
    return reg_wave_data[0].freq_list[reg_wave_data[0].freq_list.size() - 1];
}

/*******************************************************************************
 * HydroData::GetNumFreqs()
 * returns number of frequencies computed
 *******************************************************************************/
double HydroData::GetNumFreqs() const {
    return reg_wave_data[0].freq_list.size();
}

/*******************************************************************************
 * HydroData::GetExcitationMagInterp()
 * returns excitation magnitudes for body b, row i, column j, frequency ix k
 *******************************************************************************/
double HydroData::GetExcitationMagInterp(int b, int i, int j, double freq_index_des) const {
    double freq_interp_val    = freq_index_des - floor(freq_index_des);
    double excitationMagFloor = GetExcitationMagVal(b, i, j, floor(freq_index_des));
    double excitationMagCeil  = GetExcitationMagVal(b, i, j, floor(freq_index_des) + 1);
    double excitationMag      = (freq_interp_val * (excitationMagCeil - excitationMagFloor)) + excitationMagFloor;

    return excitationMag;
}

/*******************************************************************************
 * HydroData::GetExcitationMagValue()
 * returns excitation magnitudes for row i, column j, frequency ix k
 *******************************************************************************/
double HydroData::GetExcitationMagVal(int b, int i, int j, int k) const {
    int indexExMag = k + reg_wave_data[b].excitation_mag_matrix.dimension(2) * i;
    return reg_wave_data[b].excitation_mag_matrix(i, j, k) * sim_data.rho * sim_data.g;
}

/*******************************************************************************
 * HydroData::GetExcitationPhaseValue()
 * returns excitation phases for row i, column j, frequency k
 *******************************************************************************/
double HydroData::GetExcitationPhaseVal(int b, int i, int j, int k) const {
    return reg_wave_data[b].excitation_phase_matrix(i, j, k);
}

/*******************************************************************************
 * HydroData::GetExcitationPhaseInterp()
 * returns excitation phases for row i, column j, frequency ix k
 *******************************************************************************/
double HydroData::GetExcitationPhaseInterp(int b, int i, int j, double freq_index_des) const {
    double freq_interp_val      = freq_index_des - floor(freq_index_des);
    double excitationPhaseFloor = GetExcitationPhaseVal(b, i, j, floor(freq_index_des));
    double excitationPhaseCeil  = GetExcitationPhaseVal(b, i, j, floor(freq_index_des) + 1);
    double excitationPhase = (freq_interp_val * (excitationPhaseCeil - excitationPhaseFloor)) + excitationPhaseFloor;

    return excitationPhase;
}

/*******************************************************************************
 * H5FileInfo::GetExcitationIRFval()
 * returns rirf val for DoF: 0,...,5; col: 0,...,6N-1; s: 0,...,1001 excitation_irf_dims[2]
 *******************************************************************************/
double H5FileInfo::GetExcitationIRFVal(int dof, int col, int s) const {
    int index = s + excitation_irf_dims[2] * (col + dof * excitation_irf_dims[1]);  // TODO check index
    if (index < 0 || index >= excitation_irf_dims[0] * rirf_dims[1] * excitation_irf_dims[2]) {
        std::cout << "out of bounds IRF\n";
        return 0;
    } else {
        return excitation_irf_matrix[index] * _rho * _g;  // scale radiation force by rho
    }
}

/*******************************************************************************
 * H5FileInfo::GetExcitationIRFTime()
 * returns the std::vector of rirf_time_vector from h5 file
 *******************************************************************************/
std::vector<double> H5FileInfo::GetExcitationIRFTime() const {
    return excitation_irf_time;
}

/*******************************************************************************
 * H5FileInfo::GetExcitationIRF()
 * returns the std::vector of excitation_irf_matrix from h5 file
 *******************************************************************************/
std::vector<double> H5FileInfo::GetExcitationIRF() const {
    return excitation_irf_matrix;
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> ResampleTimeSeries(const Eigen::VectorXd& time_series,
                                                               double dt_old,
                                                               double dt_new) {
    if (dt_new == dt_old) {
        // If the new time resolution is the same as the original, return the original time series
        Eigen::VectorXd t_old = Eigen::VectorXd::LinSpaced(time_series.size(), 0, (time_series.size() - 1) * dt_old);
        t_old.array() -= 0.5 * t_old[t_old.size() - 1];
        return {t_old, time_series};
    }

    // Calculate the new time vector
    int newSize           = static_cast<int>(ceil(time_series.size() * dt_old / dt_new));
    Eigen::VectorXd t_new = Eigen::VectorXd::LinSpaced(newSize, 0, (time_series.size() - 1) * dt_old);

    // Create the original time vector
    Eigen::VectorXd t_old = Eigen::VectorXd::LinSpaced(time_series.size(), 0, (time_series.size() - 1) * dt_old);

    Eigen::VectorXd time_series_new(newSize);

    // Interpolate using cubic spline interpolation
    Eigen::Spline<double, 1> spline =
        Eigen::SplineFitting<Eigen::Spline<double, 1>>::Interpolate(time_series.transpose(), 3, t_old);

    for (int i = 0; i < newSize; i++) {
        time_series_new[i] = spline(t_new[i])[0];
    }

    // Shift t_new to the left by half of the max value of t_old
    t_new.array() -= 0.5 * t_old[t_old.size() - 1];

    return {t_new, time_series_new};
}

void H5FileInfo::ResampleExcitationIRFTime(double dt_new) {
    Eigen::VectorXd excitation_irf_t(excitation_irf_time.size());
    for (size_t i = 0; i < excitation_irf_time.size(); i++) {
        excitation_irf_t[i] = excitation_irf_time[i];
    }
    double excitation_irf_dt = excitation_irf_time[1] - excitation_irf_time[0];
    std::pair<Eigen::VectorXd, Eigen::VectorXd> resampled_excitation_irf_time =
        ResampleTimeSeries(excitation_irf_t, excitation_irf_dt, dt_new);
    excitation_irf_time_resampled    = resampled_excitation_irf_time.first;
    is_excitation_irf_time_resampled = true;
}

// std::pair<Eigen::VectorXd, Eigen::VectorXd> H5FileInfo::ResampleExcitationIRF(double dt_new) {
//    Eigen::VectorXd excitation_irf(excitation_irf_matrix.size());
//    for (size_t i = 0; i < excitation_irf_matrix.size(); i++) {
//        excitation_irf[i] = excitation_irf_matrix[i];
//    }
//    double excitation_irf_dt = excitation_irf_time[1] - excitation_irf_time[0];
//    return ResampleTimeSeries(excitation_irf, excitation_irf_dt, dt_new);
//}

void H5FileInfo::ResampleExcitationIRF(double dt_new) {
    Eigen::VectorXd excitation_irf(excitation_irf_matrix.size());
    for (size_t i = 0; i < excitation_irf_matrix.size(); i++) {
        excitation_irf[i] = excitation_irf_matrix[i];
    }
    double excitation_irf_dt = excitation_irf_time[1] - excitation_irf_time[0];
    std::pair<Eigen::VectorXd, Eigen::VectorXd> resampled_excitation_irf =
        ResampleTimeSeries(excitation_irf, excitation_irf_dt, dt_new);
    excitation_irf_time_resampled    = resampled_excitation_irf.first;
    excitation_irf_resampled         = resampled_excitation_irf.second;
    is_excitation_irf_time_resampled = true;
    is_excitation_irf_resampled      = true;
}

/*******************************************************************************
 * H5FileInfo::GetExcitationIRFTimeResampled()
 * returns the Eigen::VectorXd of excitation_irf_time_resampled (ResampleExcitationIRFTime() needs to be called first)
 *******************************************************************************/
Eigen::VectorXd H5FileInfo::GetExcitationIRFTimeResampled() const {
    if (!is_excitation_irf_time_resampled) {
        std::cerr << "Warning: ResampleExcitationIRFTime() has not been called before accessing resampled data. "
                     "Returning original excitation IRF time vector.\n";
        Eigen::VectorXd excitation_irf_time_eigen(excitation_irf_time.size());
        for (size_t i = 0; i < excitation_irf_time.size(); i++) {
            excitation_irf_time_eigen[i] = excitation_irf_time[i];
        }
        return excitation_irf_time_eigen;
    }
    return excitation_irf_time_resampled;
}

/*******************************************************************************
 * H5FileInfo::GetExcitationIRFResampled()
 * returns the Eigen::VectorXd of excitation_irf_time_resampled (ResampleExcitationIRFTime() needs to be called first)
 *******************************************************************************/
Eigen::VectorXd H5FileInfo::GetExcitationIRFResampled() const {
    if (!is_excitation_irf_resampled) {
        std::cerr << "Warning: ResampleExcitationIRFTime() has not been called before accessing resampled data. "
                     "Returning original excitation IRF time vector.\n";
        Eigen::VectorXd excitation_irf_eigen(excitation_irf_matrix.size());
        for (size_t i = 0; i < excitation_irf_matrix.size(); i++) {
            excitation_irf_eigen[i] = excitation_irf_matrix[i];
        }
        return excitation_irf_eigen;
    }
    return excitation_irf_resampled;
}

/*******************************************************************************
 * H5FileInfo::GetExcitationIRFResampledVal()
 * returns rirf val for DoF: 0,...,5; col: 0,...,6N-1; s: 0,...,1001 excitation_irf_dims[2]
 *******************************************************************************/
double H5FileInfo::GetExcitationIRFResampledVal(int dof, int col, int s) const {
    int index = s + excitation_irf_dims[2] * (col + dof * excitation_irf_dims[1]);  // TODO check index
    if (index < 0 || index >= excitation_irf_dims[0] * rirf_dims[1] * excitation_irf_dims[2]) {
        std::cout << "out of bounds IRF\n";
        return 0;
    } else {
        return excitation_irf_matrix[index] * _rho * _g;  // scale radiation force by rho
    }
}