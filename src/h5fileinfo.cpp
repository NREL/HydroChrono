#include <hydroc/h5fileinfo.h>

#include <filesystem>

#include <H5Cpp.h>

using namespace chrono;

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
H5FileInfo::H5FileInfo(std::string file, std::string Name) {
	h5_file_name = file;
	bodyName = Name;
	std::cout << "searching for file: " << file << std::endl;
	if (std::filesystem::exists(file)) {
		std::cout << "found file at: " << std::filesystem::absolute(file) << std::endl;
	}
	else {
		std::cout << "h5 file does not exist, absolute file location: " << std::filesystem::absolute(file) << std::endl;
	}
	readH5Data();
}

/*******************************************************************************
* H5FileInfo copy constructor (H5FileInfo& rhs)
* defines basic copy constructor using the = operator
*******************************************************************************/
H5FileInfo::H5FileInfo(H5FileInfo& old) {
	*this = old;
}

/*******************************************************************************
* H5FileInfo::operator = (H5FileInfo& rhs) 
* defines = operator, sets the left object = to right object (*this = rhs)
* returns *this
*******************************************************************************/
H5FileInfo& H5FileInfo::operator = (H5FileInfo& rhs) {
	printed = rhs.printed;
	cg = rhs.cg;
	cb = rhs.cb;
	bodyNum = rhs.bodyNum;
	_rho = rhs._rho;
	_g = rhs._g;
	_disp_vol = rhs._disp_vol;
	freq_list = rhs.freq_list;
	lin_matrix = rhs.lin_matrix;
	inf_added_mass = rhs.inf_added_mass;
	rirf_matrix = rhs.rirf_matrix;
	rirf_dims = rhs.rirf_dims;
	rirf_time_vector = rhs.rirf_time_vector;
	radiation_damping_matrix = rhs.radiation_damping_matrix;
	Bw_dims = rhs.Bw_dims;
	excitation_mag_matrix = rhs.excitation_mag_matrix;
	excitation_mag_dims = rhs.excitation_mag_dims;
	excitation_re_matrix = rhs.excitation_re_matrix;
	re_dims = rhs.re_dims;
	excitation_im_matrix = rhs.excitation_im_matrix;
	im_dims = rhs.im_dims;
	excitation_phase_matrix = rhs.excitation_phase_matrix;
	excitation_phase_dims = rhs.excitation_phase_dims;
	h5_file_name = rhs.h5_file_name;
	bodyName = rhs.bodyName;
	return *this;
}

/*******************************************************************************
* H5FileInfo::readH5Data()
* private member function called from constructor
* calls Initialize functions to read h5 file information into  member variables
*******************************************************************************/
void H5FileInfo::readH5Data() {
	// open file with read only access
	H5::H5File userH5File(h5_file_name, H5F_ACC_RDONLY);

	InitScalar(userH5File, "simulation_parameters/rho", _rho);
	InitScalar(userH5File, "simulation_parameters/g", _g);
	InitScalar(userH5File, bodyName + "/properties/disp_vol", _disp_vol);

	Init1D(userH5File, bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/t", rirf_time_vector);
	Init1D(userH5File, bodyName + "/properties/cb", cb);
	Init1D(userH5File, bodyName + "/properties/cg", cg);
	Init1D(userH5File, "simulation_parameters/w", freq_list);

	Init2D(userH5File, bodyName + "/hydro_coeffs/linear_restoring_stiffness", lin_matrix);
	Init2D(userH5File, bodyName + "/hydro_coeffs/added_mass/inf_freq", inf_added_mass);

	Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/mag", excitation_mag_matrix, excitation_mag_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/re", excitation_re_matrix, re_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/im", excitation_im_matrix, im_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/K", rirf_matrix, rirf_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/radiation_damping/all", radiation_damping_matrix, Bw_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/phase", excitation_phase_matrix, excitation_phase_dims);
	
	// use same scalar function to set the int valued body number
	double temp;
	InitScalar(userH5File, bodyName + "/properties/body_number", temp);
	bodyNum = (int)temp;

	//_rirf_timestep = rirf_time_vector[1] - rirf_time_vector[0]; //N.B. assumes RIRF has fixed timestep.
	userH5File.close();
}

/*******************************************************************************
* H5FileInfo::InitScalar(file, data_name, var)
* file - H5File& to open h5 file with body information
* data_name - string name of H5::DataSet to be read within h5 file
* var - member variable to store the information in data_name in
* Reads a double type variable from data_name DataSet in file, stores it in var
*******************************************************************************/
void H5FileInfo::InitScalar(H5::H5File& file, std::string data_name, double& var) {
	H5::DataSet dataset = file.openDataSet(data_name);
	H5::DataSpace filespace = dataset.getSpace();
	hsize_t dims[2] = { 0,0 };
	int rank = filespace.getSimpleExtentDims(dims);
	H5::DataSpace mspace1 = H5::DataSpace(rank, dims);
	double* temp = new double[dims[0] * dims[1]];
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
void H5FileInfo::Init1D(H5::H5File& file, std::string data_name, std::vector<double>& var) {
	// open specific dataset
	H5::DataSet dataset = file.openDataSet(data_name);
	// Get filespace for rank and dimension
	H5::DataSpace filespace = dataset.getSpace();
	// Get number of dimensions in the file dataspace
	// Get and print the dimension sizes of the file dataspace
	hsize_t dims[2] = { 0,0 };    // dataset dimensions
	int rank = filespace.getSimpleExtentDims(dims);
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
void H5FileInfo::Init2D(H5::H5File& file, std::string data_name, ChMatrixDynamic<double>& var) {
	//data_name = bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/K";
	H5::DataSet dataset = file.openDataSet(data_name);
	H5::DataSpace filespace = dataset.getSpace();
	hsize_t dims[2] = { 0,0 };
	int rank = filespace.getSimpleExtentDims(dims);
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
void H5FileInfo::Init3D(H5::H5File& file, std::string data_name, std::vector<double>& var, std::vector<int>& d) {
	// open specific dataset
	H5::DataSet dataset = file.openDataSet(data_name);
	// Get filespace for rank and dimension
	H5::DataSpace filespace = dataset.getSpace();
	hsize_t dims[3] = { 0,0,0 };
	int rank = filespace.getSimpleExtentDims(dims);
	// read file into data_out 2d array
	H5::DataSpace mspace(rank, dims);
	// rirf_dims[0] is number of rows, rirf_dims[1] is number of columns, rirf_dims[2] is number of matrices
	double* temp = new double[dims[0] * dims[1] * dims[2]];
	// read file info into data_out, a 2d array
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace, filespace);
	// set var here
	var.resize(dims[0] * dims[1] * dims[2]);
	d.resize(3);
	for (int i = 0; i < 3; i++) {
		d[i] = dims[i];
	}
	for (int i = 0; i < dims[0] * dims[1] * dims[2]; i++) {
		var[i] = temp[i];
	}
	dataset.close();
	delete[] temp;
}

/*******************************************************************************
* H5FileInfo::~H5FileInfo()
* H5FileInfo destructor TODO
*******************************************************************************/
H5FileInfo::~H5FileInfo() { }

/*******************************************************************************
* H5FileInfo::GetRIRFDims(int i) returns the i-th component of the dimensions of radiation_damping_matrix
* i = [0,1,2] -> [number of rows, number of columns, number of matrices]
*******************************************************************************/
int H5FileInfo::GetRIRFDims(int i) const { 
	return rirf_dims[i];
}

/*******************************************************************************
* H5FileInfo::GetHydrostaticStiffness()
* returns the linear restoring stiffness matrix element in row i , column j
*******************************************************************************/
double H5FileInfo::GetHydrostaticStiffness(int i, int j) const { 
	return lin_matrix(i, j) * rho * g;
}


/*******************************************************************************
* H5FileInfo::GetRIRFval()
* returns rirf val for DoF: 0,...,5; col: 0,...,6N-1; s: 0,...,1001 rirfdims[2]
*******************************************************************************/
double H5FileInfo::GetRIRFval(int dof, int col, int s) const {
	int index = s + rirf_dims[2] * (col + dof * rirf_dims[1]); // TODO check index
	if (index < 0 || index >= rirf_dims[0] * rirf_dims[1] * rirf_dims[2]) {
		std::cout << "out of bounds IRF\n";
		return 0;
	}
	else {
		return rirf_matrix[index] * rho; // scale radiation force by rho
	}
}

/*******************************************************************************
* H5FileInfo::GetRIRFTimeVector()
* returns the std::vector of rirf_time_vector from h5 file
*******************************************************************************/
std::vector<double> H5FileInfo::GetRIRFTimeVector() const { 
	return rirf_time_vector;
}

/*******************************************************************************
* H5FileInfo::GetInfAddedMassMatrix()
* returns the matrix for added mass at infinite frequency scaled by rho
*******************************************************************************/
ChMatrixDynamic<double> H5FileInfo::GetInfAddedMassMatrix() const { 
	return inf_added_mass * rho;
}


/*******************************************************************************
* H5FileInfo::GetNumFreqs()
* returns number of frequencies computed
*******************************************************************************/
double H5FileInfo::GetNumFreqs() const { 
	return freq_list.size();
}

/*******************************************************************************
* H5FileInfo::GetOmegaMin()
* returns min value of omega
*******************************************************************************/
//double H5FileInfo::GetOmegaMin() const { //TODO cut this func????
//	return freq_list[0];
//}

/*******************************************************************************
* H5FileInfo::GetOmegaMax()
* returns max value of omega
*******************************************************************************/
double H5FileInfo::GetOmegaMax() const {
	return freq_list[freq_list.size() - 1];
}

/*******************************************************************************
* H5FileInfo::GetOmegaDelta()
* returns omega step size
*******************************************************************************/
double H5FileInfo::GetOmegaDelta() const {
	return GetOmegaMax() / GetNumFreqs();
}

/*******************************************************************************
* H5FileInfo::GetExcitationMagValue()
* returns excitation magnitudes for row i, column j, frequency ix k
*******************************************************************************/
double H5FileInfo::GetExcitationMagValue(int i, int j, int k) const {
	int indexExMag = k + excitation_mag_dims[2] * i;
	return excitation_mag_matrix[indexExMag] * rho * g;
}

/*******************************************************************************
* H5FileInfo::GetExcitationMagInterp()
* returns excitation magnitudes for row i, column j, frequency ix k
*******************************************************************************/
double H5FileInfo::GetExcitationMagInterp(int i, int j, double freq_index_des) const {
	double freq_interp_val = freq_index_des - floor(freq_index_des);
	double excitationMagFloor = GetExcitationMagValue(i, j, floor(freq_index_des));
	double excitationMagCeil = GetExcitationMagValue(i, j, floor(freq_index_des) +1);
	double excitationMag = (freq_interp_val * (excitationMagCeil - excitationMagFloor)) + excitationMagFloor;

	return excitationMag;
}

/*******************************************************************************
* H5FileInfo::GetExcitationPhaseValue()
* returns excitation phases for row i, column j, frequency k
*******************************************************************************/
double H5FileInfo::GetExcitationPhaseValue(int i, int j, int k) const {
	int indexExPhase = k + excitation_phase_dims[2] * i;
	return excitation_phase_matrix[indexExPhase];
}

/*******************************************************************************
* H5FileInfo::GetExcitationPhaseInterp()
* returns excitation phases for row i, column j, frequency ix k
*******************************************************************************/
double H5FileInfo::GetExcitationPhaseInterp(int i, int j, double freq_index_des) const {
	double freq_interp_val = freq_index_des - floor(freq_index_des);
	double excitationPhaseFloor = GetExcitationPhaseValue(i, j, floor(freq_index_des));
	double excitationPhaseCeil = GetExcitationPhaseValue(i, j, floor(freq_index_des) + 1);
	double excitationPhase = (freq_interp_val * (excitationPhaseCeil - excitationPhaseFloor)) + excitationPhaseFloor;

	return excitationPhase;
}
