#include "hydro_forces.h"

// =============================================================================
// H5FileInfo Class Definitions
// =============================================================================

/*******************************************************************************
* H5FileInfo::read_data()
* private member function called from constructor
* reads h5 file data and stores it in member variables for use with other
* classes and forces
*******************************************************************************/
void H5FileInfo::read_data() {
	// testing HDF5 compatibility
	// open file with read only access
	H5::H5File sphereFile(h5_file_name, H5F_ACC_RDONLY); // "../../test_for_chrono/sphere.h5"

	//
	// Read linear restoring stiffness file info into matrices
	// "body1"
	std::string data_name = bodyNum + "/hydro_coeffs/linear_restoring_stiffness"; // "body1/hydro_coeffs/linear_restoring_stiffness"
	H5::DataSet dataset = sphereFile.openDataSet(data_name);
	// Get filespace for rank and dimension
	H5::DataSpace filespace = dataset.getSpace();
	// Get number of dimensions in the file dataspace
	// Get and print the dimension sizes of the file dataspace
	hsize_t dims[2];    // dataset dimensions
	int rank = filespace.getSimpleExtentDims(dims);
	// read file into data_out 2d array
	H5::DataSpace mspace1(rank, dims);
	double *temp;
	temp = new double[dims[0] * dims[1]];
	// read file info into data_out, a 2d array
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// turn the 2d array into a ChMatrix (Eigen dynamic matrix)
	lin_matrix.resize(dims[0], dims[1]);
	for (int i = 0; i < dims[0]; i++) {
		for (int j = 0; j < dims[1]; j++) {
			lin_matrix(i, j) = temp[i * dims[1] + j];
		}
	}
	dataset.close();
	delete[] temp;

	data_name = bodyNum + "/hydro_coeffs/added_mass/inf_freq";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	inf_added_mass.resize(dims[0], dims[1]);
	for (int i = 0; i < dims[0]; i++) {
		for (int j = 0; j < dims[1]; j++) {
			inf_added_mass(i, j) = temp[i * dims[1] + j];
		}
	}
	dataset.close();
	delete[] temp;

	// repeat same steps from above to get the cb and cg...reusing some of the previous arrays etc
	data_name = bodyNum + "/properties/cb";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	for (int i = 0; i < dims[0]; i++) {
		cb[i] = temp[i];
	}
	dataset.close();

	// repeat finally for cg
	data_name = bodyNum + "/properties/cg";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	for (int i = 0; i < dims[0]; i++) {
		cg[i] = temp[i];
	}
	dataset.close();
	delete[] temp;

	// read displaced volume for buoyancy force
	data_name = bodyNum + "/properties/disp_vol";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	disp_vol = temp[0];
	dataset.close();

	// read rho
	data_name = "simulation_parameters/rho";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	rho = temp[0];
	dataset.close();

	// read g
	data_name = "simulation_parameters/g";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	g = temp[0];
	dataset.close();
	delete[] temp;
	lin_matrix *= rho*g; // scale by rho*g

	// read frequencies
	data_name = "simulation_parameters/w";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(freq_dims);
	mspace1 = H5::DataSpace(rank, freq_dims);
	temp = new double[freq_dims[0] * freq_dims[1]];
	freq_list.resize(freq_dims[0]);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into timesteps chvector
	for (int i = 0; i < freq_dims[0]; i++) {
		freq_list[i] = temp[i];
	}
	dataset.close();
	delete[] temp;

	// read K
	data_name = bodyNum + "/hydro_coeffs/radiation_damping/impulse_response_fun/K";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	int rank3 = filespace.getSimpleExtentDims(rirf_dims);
	// read file into data_out 2d array
	H5::DataSpace mspace3(rank3, rirf_dims);
	// rirf_dims[0] is number of rows, rirf_dims[1] is number of columns, rirf_dims[2] is number of matrices
	rirf_matrix = new double[rirf_dims[0] * rirf_dims[1] * rirf_dims[2]];
	// read file info into data_out, a 2d array
	dataset.read(rirf_matrix, H5::PredType::NATIVE_DOUBLE, mspace3, filespace);
	dataset.close();

	// read timesteps vector
	data_name = bodyNum + "/hydro_coeffs/radiation_damping/impulse_response_fun/t";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	timesteps.resize(dims[0]);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into timesteps chvector
	for (int i = 0; i < dims[0]; i++) {
		timesteps[i] = temp[i];
	}
	dataset.close();
	delete[] temp;

	// read B(w)
	data_name = bodyNum + "/hydro_coeffs/radiation_damping/all";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank3 = filespace.getSimpleExtentDims(radiation_damping_dims);
	// read file into data_out 2d array
	H5::DataSpace mspaceB(rank3, radiation_damping_dims);
	// radiation_damping_dims[0] is number of rows, radiation_damping_dims[1] is number of columns, radiation_damping_dims[2] is number of matrices
	radiation_damping_matrix = new double[radiation_damping_dims[0] * radiation_damping_dims[1] * radiation_damping_dims[2]];
	// read file info into data_out, a 2d array
	dataset.read(radiation_damping_matrix, H5::PredType::NATIVE_DOUBLE, mspaceB, filespace);
	dataset.close();

	// read excitation force coefficients - magnitude
	data_name = bodyNum + "/hydro_coeffs/excitation/mag";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank3 = filespace.getSimpleExtentDims(excitation_mag_dims);
	// read file into data_out 2d array
	H5::DataSpace mspaceExMag(rank3, excitation_mag_dims);
	//// excitation_mag_dims[0] is number of rows, excitation_mag_dims[1] is number of columns, excitation_mag_dims[2] is number of matrices
	excitation_mag_matrix = new double[excitation_mag_dims[0] * excitation_mag_dims[1] * excitation_mag_dims[2]];
	//// read file info into data_out, a 2d array
	dataset.read(excitation_mag_matrix, H5::PredType::NATIVE_DOUBLE, mspaceExMag, filespace);
	dataset.close();

	// read excitation force coefficients - phase
	data_name = bodyNum + "/hydro_coeffs/excitation/phase";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank3 = filespace.getSimpleExtentDims(excitation_phase_dims);
	// read file into data_out 2d array
	H5::DataSpace mspaceExPhase(rank3, excitation_phase_dims);
	//// excitation_phase_dims[0] is number of rows, excitation_phase_dims[1] is number of columns, excitation_phase_dims[2] is number of matrices
	excitation_phase_matrix = new double[excitation_phase_dims[0] * excitation_phase_dims[1] * excitation_phase_dims[2]];
	//// read file info into data_out, a 2d array
	dataset.read(excitation_phase_matrix, H5::PredType::NATIVE_DOUBLE, mspaceExPhase, filespace);
	dataset.close();

	// read excitation force coefficients - real
	data_name = bodyNum + "/hydro_coeffs/excitation/re";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank3 = filespace.getSimpleExtentDims(excitation_re_dims);
	// read file into data_out 2d array
	H5::DataSpace mspaceExRe(rank3, excitation_re_dims);
	//// excitation_phase_dims[0] is number of rows, excitation_phase_dims[1] is number of columns, excitation_phase_dims[2] is number of matrices
	excitation_re_matrix = new double[excitation_re_dims[0] * excitation_re_dims[1] * excitation_re_dims[2]];
	//// read file info into data_out, a 2d array
	dataset.read(excitation_re_matrix, H5::PredType::NATIVE_DOUBLE, mspaceExRe, filespace);
	dataset.close();

	// read excitation force coefficients - im
	data_name = bodyNum + "/hydro_coeffs/excitation/im";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank3 = filespace.getSimpleExtentDims(excitation_im_dims);
	// read file into data_out 2d array
	H5::DataSpace mspaceExIm(rank3, excitation_im_dims);
	//// excitation_phase_dims[0] is number of rows, excitation_phase_dims[1] is number of columns, excitation_phase_dims[2] is number of matrices
	excitation_im_matrix = new double[excitation_im_dims[0] * excitation_im_dims[1] * excitation_im_dims[2]];
	//// read file info into data_out, a 2d array
	dataset.read(excitation_im_matrix, H5::PredType::NATIVE_DOUBLE, mspaceExIm, filespace);
	dataset.close();

	sphereFile.close();
}

H5FileInfo::H5FileInfo() {}

H5FileInfo::~H5FileInfo() {
	delete[] rirf_matrix;
}

/*******************************************************************************
* H5FileInfo constructor
* requires file name (in absolute file name or referenced from executable location)
* and body name (name of body's section in H5 file, ie "body1" etc)
* each body in system should have its own H5FileInfo object
* calls read_data()
*******************************************************************************/
H5FileInfo::H5FileInfo(std::string file, std::string bodyName) {
	h5_file_name = file;
	bodyNum = bodyName;
	read_data();
}

/*******************************************************************************
* H5FileInfo::get_lin_matrix()
* returns the linear restoring stiffness matrix
*******************************************************************************/
ChMatrixDynamic<double> H5FileInfo::get_lin_matrix() const {
	return lin_matrix;
}

/*******************************************************************************
* H5FileInfo::get_inf_added_mass_matrix()
* returns the added mass matrix at infinite frequency
*******************************************************************************/
ChMatrixDynamic<double> H5FileInfo::get_inf_added_mass_matrix() const {
	return inf_added_mass * rho;
}

/*******************************************************************************
* H5FileInfo::get_equil_cg()
* returns cg, center of gravity of object's body
*******************************************************************************/
ChVector<> H5FileInfo::get_equil_cg() const {
	return cg;
}

/*******************************************************************************
* H5FileInfo::get_equil_cb()
* returns cb, the rotational equilibrium
*******************************************************************************/
ChVector<> H5FileInfo::get_equil_cb() const {
	return cb;
}

/*******************************************************************************
* H5FileInfo::get_rho()
* returns the density of water, rho (kg/m^3 usually)
*******************************************************************************/
double H5FileInfo::get_rho() const {
	return rho;
}

/*******************************************************************************
* H5FileInfo::get_g()
* returns g, gravitational acceleration, m/s^2
*******************************************************************************/
double H5FileInfo::get_g() const {
	return g;
}

/*******************************************************************************
* H5FileInfo::get_disp_vol()
* returns displaced volume when body at equilibrium, m^3
*******************************************************************************/
double H5FileInfo::get_disp_vol() const {
	return disp_vol;
}

/*******************************************************************************
* H5FileInfo::get_impulse_resp_matrix()
* returns impulse response coeff for row m, column n, step s
*******************************************************************************/
double H5FileInfo::get_rirf_val(int m, int n, int s) const {
	int index = s + rirf_dims[2] * (n + m * rirf_dims[1]);
	if (index < 0 || index >= rirf_dims[0] * rirf_dims[1] * rirf_dims[2]) {
		std::cout << "out of bounds IRF\n";
		return 0;
	}
	else {
		return rirf_matrix[index] * get_rho(); // scale radiation force by rho
	}
}

/*******************************************************************************
* H5FileInfo::get_rirf_dims(int i) returns the i-th component of the dimensions of radiation_damping_matrix
* i = [0,1,2] -> [number of rows, number of columns, number of matrices]
*******************************************************************************/
int H5FileInfo::get_rirf_dims(int i) const {
	return rirf_dims[i];
}

//TODO: Get B(w)

/*******************************************************************************
* H5FileInfo::get_num_freqs()
* returns number of frequencies computed
*******************************************************************************/
double H5FileInfo::get_num_freqs() const {
	return freq_dims[0];
}

/*******************************************************************************
* H5FileInfo::get_omega_min()
* returns min value of omega 
*******************************************************************************/
double H5FileInfo::get_omega_min() const {
	return freq_list[0];
}

/*******************************************************************************
* H5FileInfo::get_omega_max()
* returns max value of omega
*******************************************************************************/
double H5FileInfo::get_omega_max() const {
	return freq_list[freq_dims[0] - 1];
}

/*******************************************************************************
* H5FileInfo::get_domega()
* returns omega step size
*******************************************************************************/
double H5FileInfo::get_domega() const {
	return get_omega_max() / get_num_freqs();
}

/*******************************************************************************
* H5FileInfo::get_excitation_mag_ix()
* returns excitation magnitudes for row i, column j, frequency ix k
*******************************************************************************/
double H5FileInfo::get_excitation_mag_ix(int i, int j, int k) const {
	int indexExMag = k + excitation_mag_dims[2] * i;
	return excitation_mag_matrix[indexExMag] * rho * g;
}

/*******************************************************************************
* H5FileInfo::get_excitation_mag_interp()
* returns excitation magnitudes for row i, column j, frequency ix k
*******************************************************************************/
double H5FileInfo::get_excitation_mag_interp(int i, int j, double freqIndexDes) const {
	double freqInterpVal = freqIndexDes - floor(freqIndexDes);
	double excitationMagFloor = get_excitation_mag_ix(i, j, floor(freqIndexDes));
	double excitationMagCeil = get_excitation_mag_ix(i, j, floor(freqIndexDes) +1);
	double excitationMag = (freqInterpVal * (excitationMagCeil - excitationMagFloor)) + excitationMagFloor;

	return excitationMag;
}

/*******************************************************************************
* H5FileInfo::get_excitation_phase_ix()
* returns excitation phases for row i, column j, frequency k
*******************************************************************************/
double H5FileInfo::get_excitation_phase_ix(int i, int j, int k) const {
	int indexExPhase = k + excitation_phase_dims[2] * i;
	return excitation_phase_matrix[indexExPhase];
}

/*******************************************************************************
* H5FileInfo::get_excitation_phase_interp()
* returns excitation phases for row i, column j, frequency ix k
*******************************************************************************/
double H5FileInfo::get_excitation_phase_interp(int i, int j, double freqIndexDes) const {
	double freqInterpVal = freqIndexDes - floor(freqIndexDes);
	double excitationPhaseFloor = get_excitation_phase_ix(i, j, floor(freqIndexDes));
	double excitationPhaseCeil = get_excitation_phase_ix(i, j, floor(freqIndexDes) + 1);
	double excitationPhase = (freqInterpVal * (excitationPhaseCeil - excitationPhaseFloor)) + excitationPhaseFloor;

	return excitationPhase;
}

/*******************************************************************************
* H5FileInfo::get_delta_t() returns the difference in first 2 timesteps
*******************************************************************************/
double H5FileInfo::get_delta_t() const {
	return timesteps[1] - timesteps[0];
}

/*******************************************************************************
* H5FileInfo::get_times()
* returns the vector of timesteps from h5 file
*******************************************************************************/
std::vector<double> H5FileInfo::get_times() const {
	return timesteps;
}

// =============================================================================
// ForceTorqueFunc Class Definitions
// =============================================================================

/*******************************************************************************
* ForceTorqueFunc constructor
* sets pointer to HydroForces member object and index for which component
* this ForceTorqueFunc object refers to
*******************************************************************************/
ForceTorqueFunc::ForceTorqueFunc(HydroForces* b, int i) : base(b), index(i) {	}

/*******************************************************************************
* ForceTorqueFunc::Clone()
* required override function since ForceTorqueFunc inherits from ChFunction
*******************************************************************************/
ForceTorqueFunc* ForceTorqueFunc::Clone() const  {
	return new ForceTorqueFunc(*this);
}

/*******************************************************************************
* ForceTorqueFunc::Get_y()
* returns the value of the index-th coordinate of the linear restoring force
* at each time
* required override function since ForceTorqueFunc inherits from ChFunction
*******************************************************************************/
double ForceTorqueFunc::Get_y(double x) const {
	return base->coordinateFunc(index);
}

/*******************************************************************************
* ForceTorqueFunc::SetBase
* currently unused(?) function to set base if not set in constructor
*******************************************************************************/
void ForceTorqueFunc::SetBase(HydroForces* b) {
	base = b;
}

/*******************************************************************************
* ForceTorqueFunc::SetIndex
* currently unused(?) function to set index if not set in constructor
*******************************************************************************/
void ForceTorqueFunc::SetIndex(int i) {
	index = i;
}

// =============================================================================
// HydroInputs Class Definitions
// =============================================================================

HydroInputs::HydroInputs() {
	
}

	
//public:
//	HydroInputs();
//	~HydroInputs();
//	double regularWaveAmplitude;
//	double regularWavePeriod;
//	double regularWaveOmega;
//	double get_regular_wave_omega(double regularWavePeriod);
//private:
//	
//};

// =============================================================================
// HydroForces Class Definitions
// =============================================================================

/*******************************************************************************
* HydroForces default constructor
* initializes array of ForceTorqueFunc objects and pointers to each force/torque
*******************************************************************************/
HydroForces::HydroForces() : forces{ {this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5} } {
	for (unsigned i = 0; i < 6; i++) {
		force_ptrs[i] = std::shared_ptr<ForceTorqueFunc>(forces + i, [](ForceTorqueFunc*) {});
		// sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is
		// automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
		// them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the
		// default deletion logic to do nothing
		// Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
	}
	chronoForce = chrono_types::make_shared<ChForce>();
	chronoTorque = chrono_types::make_shared<ChForce>();
}

/*******************************************************************************
* HydroForces constructor
* calls default constructor and initializes hydro force info
* from H5FileInfo
* also initializes ChBody that this force will be applied to
*******************************************************************************/
HydroForces::HydroForces(H5FileInfo& sysH5FileInfo, std::shared_ptr<ChBody> object, HydroInputs userHydroInputs) : HydroForces() {
	body = object;
	fileInfo = sysH5FileInfo;
	hydroInputs = userHydroInputs;
	// define wave inputs here
	// TODO: switch depending on wave option (regular, regularCIC, irregular, noWaveCIC)
	waveAmplitude = hydroInputs.regularWaveAmplitude;
	waveOmega = hydroInputs.regularWaveOmega;
	domega = fileInfo.get_domega();
	freqIndexDes = (waveOmega / domega) - 1;
	//int freqIndexFloor = floor(freqIndexDes);
	//double freqInterpVal = freqIndexDes - freqIndexFloor;

	for (int rowEx = 0; rowEx < 6; rowEx++) {
		forceExcitationMag[rowEx] = fileInfo.get_excitation_mag_interp(rowEx, 0, freqIndexDes); //get_excitation_mag_ix(rowEx, 0, 41);
		forceExcitationPhase[rowEx] = fileInfo.get_excitation_phase_interp(rowEx, 0, freqIndexDes);
	}

	//std::ofstream myfile;
	//myfile.open("C:\\code\\chrono_hydro_dev\\test_for_chrono_build\\Release\\freqIndex.txt");
	//myfile << freqIndex << "\n";
	//myfile.close();

	equil << fileInfo.get_equil_cg().eigen(), 0, 0, 0; // set equilibrium to (cg0, cg1, cg2, 0, 0, 0)
	prevTime = -1;
	prevTimeIRF = -1;
	velocityHistory.resize(fileInfo.get_rirf_dims(2));
	timeSteps = fileInfo.get_times();
	ChVectorN<double, 6> temp;
	for (int i = 0; i < 6; i++) {
		temp[i] = 0;
		currentForce[i] = 0;
	}
	for (int i = 0; i < 1001; i++) {
		velocityHistory[i] = temp;
	}
	offset = 0;
}

/*******************************************************************************
* HydroForces::fHydrostaticStiffness()
* calculates the matrix multiplication each time step for linear restoring stiffness
* f = [linear restoring stiffness matrix] [displacement vector]
*******************************************************************************/
ChVectorN<double, 6> HydroForces::fHydrostaticStiffness() {
	if (body->GetChTime() == prevTime) {
		return currentForce;
	}
	prevTime = body->GetChTime();
	currentForce << body->GetPos().eigen(), body->GetRot().Q_to_Euler123();// body->GetRot().Q_to_Euler123().eigen();// body->GetWvel_par();

	currentForce = currentForce - equil;

	double rollLeverArm = currentForce[3];
	double pitchLeverArm = currentForce[4];
	double yawLeverArm = currentForce[5];

	currentForce = -1 * fileInfo.get_lin_matrix() * currentForce;

	double buoyancy = fileInfo.get_rho() * fileInfo.get_g() * fileInfo.get_disp_vol();
	currentForce[2] += buoyancy;
	currentForce[3] += buoyancy * rollLeverArm;
	currentForce[4] += buoyancy * pitchLeverArm;
	currentForce[5] += buoyancy * yawLeverArm;

	return currentForce;
}

/*******************************************************************************
* HydroForces::fRadDamping()
*******************************************************************************/
ChVectorN<double, 6> HydroForces::fRadDamping() {
	// since convolutionIntegral called for each DoF each timestep, we only want to
	// calculate the vector force once each timestep. Save the prevTimeIRF
	if (body->GetChTime() == prevTimeIRF) {
		return forceRadiationDamping;
	}
	prevTimeIRF = body->GetChTime();
	int size = fileInfo.get_rirf_dims(2);
	// "shift" everything left 1
	offset--;
	if (offset < -1 * size) {
		offset += size;
	}
	int numRows = 6, numCols = 6;
	double* timeseries = new double[numRows * numCols * size];
	double* tmp_s = new double[numRows * size];
	// define shortcuts for accessing 1D arrays as 3D (or 2D) arrays
#define TIMESERIES(row,col,step) timeseries[(row)*numCols*size + (col)*size + (step)]
#define TMP_S(row,step) tmp_s[(row)*size + (step)]
	// set last entry as velocity
	for (int i = 0; i < 3; i++) {
		velocityHistory[(((size + offset) % size) + size) % size][i] = body->GetPos_dt()[i];
		velocityHistory[(((size + offset) % size) + size) % size][i + 3] = body->GetWvel_par()[i];
	}
	int vi;
	//#pragma omp parallel for
	for (int row = 0; row < 6; row++) {
		forceRadiationDamping[row] = 0.0;
		double fDampingCol = 0.0;
		//#pragma omp parallel for
		for (int col = 0; col < 6; col++) {
			for (int st = 0; st < size; st++) {
				TMP_S(row, st) = 0;
				vi = (((st + offset) % size) + size) % size; // vi takes care of circshift function from matLab
				TIMESERIES(row, col, st) = fileInfo.get_rirf_val(row, col, st) * velocityHistory[vi][col];
				TMP_S(row, st) += TIMESERIES(row, col, st);
				if (st > 0) {
					forceRadiationDamping[row] -= (TMP_S(row, st - 1) + TMP_S(row, st)) / 2.0 * (timeSteps[st] - timeSteps[st - 1ull]);
				}
			}
		}
	}
	// Deallocate memory
#undef TIMESERIES
#undef TMP_S
	delete[] timeseries;
	delete[] tmp_s;
	return forceRadiationDamping;
}

/*******************************************************************************
* HydroForces::fExcitationRegularFreq()
*
*******************************************************************************/
ChVectorN<double, 6> HydroForces::fExcitationRegularFreq() {
	if (body->GetChTime() == prevTimeEx) {
		return forceExcitation;
	}
	prevTimeEx = body->GetChTime();
	for (int rowEx = 0; rowEx < 6; rowEx++) {
		if (rowEx == 2) {
			forceExcitation[rowEx] = forceExcitationMag[rowEx] * waveAmplitude * cos(waveOmega * body->GetChTime() + forceExcitationPhase[rowEx]);
		}
		else {
			forceExcitation[rowEx] = 0.0;
		}
	}
	return forceExcitation;
}

/*******************************************************************************
* HydroForces::coordinateFunc
* if index is in [0,6] the corresponding vector component of the force vector
* is returned
* otherwise a warning is printed and the force is interpreted to be 0
*******************************************************************************/
double HydroForces::coordinateFunc(int i) {
	if (i >= 0 && i < 6) {
		double forceHs = fHydrostaticStiffness()[i];
		double forceRad = fRadDamping()[i];
		double forceExcitation = fExcitationRegularFreq()[i];
		double totalForce = forceHs + forceRad + forceExcitation;
		return totalForce;
	}
	else {
		std::cout << "wrong index" << std::endl;
		return 0;
	}
}

/*******************************************************************************
* HydroForces::SetForce
* used to initialize components of force (external ChForce pointer)
*******************************************************************************/
void HydroForces::SetForce() {
	chronoForce->SetF_x(force_ptrs[0]);
	chronoForce->SetF_y(force_ptrs[1]);
	chronoForce->SetF_z(force_ptrs[2]);
	body->AddForce(chronoForce);
}

/*******************************************************************************
* HydroForces::SetTorque
* used to initialize components of torque (external ChForce pointer with TORQUE flag set)
*******************************************************************************/
void HydroForces::SetTorque() {
	chronoTorque->SetF_x(force_ptrs[3]);
	chronoTorque->SetF_y(force_ptrs[4]);
	chronoTorque->SetF_z(force_ptrs[5]);
	chronoTorque->SetMode(ChForce::ForceType::TORQUE);
	body->AddForce(chronoTorque);
}

// =============================================================================
// ChLoadAddedMass Class Definitions
// =============================================================================

/*******************************************************************************
* ChLoadAddedMass constructor
* initializes body to have load applied to and added mass matrix from h5 file object
*******************************************************************************/
ChLoadAddedMass::ChLoadAddedMass(std::shared_ptr<ChBody> body,  ///< object to apply additional inertia to
	const H5FileInfo& file) : ChLoadCustom(body) {
	inf_added_mass_J = file.get_inf_added_mass_matrix(); //TODO switch all uses of H5FileInfo object to be like this, instead of copying the object each time?

	std::ofstream myfile;
	myfile.open("C:\\code\\chrono_hydro_dev\\HydroChrono_build\\Release\\inf_added_mass_J.txt");
	myfile << inf_added_mass_J << "\n";
	myfile.close();
}

/*******************************************************************************
* ChLoadAddedMass::ComputeJacobian()
* Computes Jacobian for load, in this case just the mass matrix is initialized
* as the added mass matrix
*******************************************************************************/
void ChLoadAddedMass::ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
	ChStateDelta* state_w,  ///< state speed to evaluate jacobians
	ChMatrixRef mK,         ///< result dQ/dx
	ChMatrixRef mR,         ///< result dQ/dv
	ChMatrixRef mM          ///< result dQ/da
) {
	//set mass matrix here

	jacobians->M = inf_added_mass_J;

	ChMatrixDynamic<double> massmat = jacobians->M;

	std::ofstream myfile2;
	myfile2.open("C:\\code\\chrono_hydro_dev\\HydroChrono_build\\Release\\massmat.txt");
	myfile2 << massmat << "\n";
	myfile2.close();

	// R gyroscopic damping matrix terms (6x6)
	// 0 for added mass
	jacobians->R.setZero();

	// K inertial stiffness matrix terms (6x6)
	// 0 for added mass
	jacobians->K.setZero();
}

/*******************************************************************************
* ChLoadAddedMass::LoadIntLoadResidual_Mv()
* Computes LoadIntLoadResidual_Mv for vector w, const c, and vector R
* Note R here is vector, and is not R gyroscopic damping matrix from ComputeJacobian
*******************************************************************************/
void ChLoadAddedMass::LoadIntLoadResidual_Mv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) {
	if (!this->jacobians)
		return;

	if (!loadable->IsSubBlockActive(0))
		return;

	// R+=c*M*a
	// segment gives the chunk of vector starting at the first argument, and going for as many elements as the second argument...
	// in this case, segment gets the 3vector starting at the 0th DOF's offset (ie 0)
	//R.segment(loadable->GetSubBlockOffset(0), 3) += c * (this->mass * (a_x + chrono::Vcross(a_w, this->c_m))).eigen();
	// in this case, segment gets the 3vector starting at the 0th DOF's + 3 offset (ie 3)
	//R.segment(loadable->GetSubBlockOffset(0) + 3, 3) += c * (this->mass * chrono::Vcross(this->c_m, a_x) + this->I * a_w).eigen();
	// since R is a vector, we can probably just do R += C*M*a with no need to separate w into a_x and a_w above
	R += c * jacobians->M * w;
}

// =============================================================================
// LoadAllHydroForces Class Definitions
// =============================================================================
/*******************************************************************************
*
*******************************************************************************/
LoadAllHydroForces::LoadAllHydroForces(std::shared_ptr<ChBody> object, std::string file, std::string bodyName, HydroInputs userHydroInputs) :
	sysFileInfo(file, bodyName), hydro_force(sysFileInfo, object, userHydroInputs) {

	my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();
	my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(object, sysFileInfo);

	object->GetSystem()->Add(my_loadcontainer);
	my_loadcontainer->Add(my_loadbodyinertia);

	hydro_force.SetForce();
	hydro_force.SetTorque();
}

//LoadAllHydroForces::~LoadAllHydroForces() {
//	/*delete sysFileInfo;*/
//	delete hydro_force;
//	delete userHydroInputs;
//	/*delete my_loadcontainer;
//	delete my_loadbodyinertia;*/
//}
