#include "hydro_forces.h"

#include <algorithm>
#include <numeric> // std::accumulate

// =============================================================================
// H5FileInfo Class Definitions
// =============================================================================

/*******************************************************************************
* H5FileInfo::readH5Data()
* private member function called from constructor
* reads h5 file data and stores it in member variables for use with other
* classes and forces
*******************************************************************************/
void H5FileInfo::readH5Data() { // TODO break up this function!
	// open file with read only access
	H5::H5File sphereFile(h5_file_name, H5F_ACC_RDONLY);

	// Read linear restoring stiffness file info into matrices
	std::string data_name = bodyName + "/hydro_coeffs/linear_restoring_stiffness";
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

	data_name = bodyName + "/hydro_coeffs/added_mass/inf_freq";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equilibrium chvector
	inf_added_mass.resize(dims[0], dims[1]);
	for (int i = 0; i < dims[0]; i++) {
		for (int j = 0; j < dims[1]; j++) {
			inf_added_mass(i, j) = temp[i * dims[1] + j];
		}
	}
	dataset.close();
	delete[] temp;

	// repeat same steps from above to get the cb and cg...reusing some of the previous arrays etc
	data_name = bodyName + "/properties/cb";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equilibrium chvector
	for (int i = 0; i < dims[0]; i++) {
		cb[i] = temp[i];
	}
	dataset.close();

	// repeat finally for cg
	// TODO add a helper function for the repeated reading things
	data_name = bodyName + "/properties/cg";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equilibrium chvector
	for (int i = 0; i < dims[0]; i++) {
		cg[i] = temp[i];
	}
	dataset.close();
	delete[] temp;

	// read displaced volume for buoyancy force
	data_name = bodyName + "/properties/disp_vol";
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

	// read body number (for index manipulation later in multibody system)
	data_name = bodyName + "/properties/body_number";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	bodyNum = (int)temp[0];
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
	// put into rirf_time_vector chvector
	for (int i = 0; i < freq_dims[0]; i++) {
		freq_list[i] = temp[i];
	}
	dataset.close();
	delete[] temp;

	// read K
	//TODO NEED new solution for this for multibody situations...
	data_name = bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/K";
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

	// read rirf_time_vector 1001x1 array
	data_name = bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/t";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	rirf_time_vector.resize(dims[0]); // resize std::vector for rirf time to be 1001
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into rirf_time_vector chvector
	for (int i = 0; i < dims[0]; i++) {
		rirf_time_vector[i] = temp[i];
	}
	rirf_timestep = rirf_time_vector[1] - rirf_time_vector[0]; //N.B. assumes RIRF has fixed timestep.
	dataset.close();
	delete[] temp;

	// read B(w)
	data_name = bodyName + "/hydro_coeffs/radiation_damping/all";
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
	data_name = bodyName + "/hydro_coeffs/excitation/mag";
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
	data_name = bodyName + "/hydro_coeffs/excitation/phase";
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
	data_name = bodyName + "/hydro_coeffs/excitation/re";
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
	data_name = bodyName + "/hydro_coeffs/excitation/im";
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
* calls readH5Data()
*******************************************************************************/
H5FileInfo::H5FileInfo(std::string file, std::string Name) {
	h5_file_name = file;
	bodyName = Name;
	readH5Data();
}

/*******************************************************************************
* H5FileInfo::GetHydrostaticStiffnessMatrix()
* returns the linear restoring stiffness matrix
*******************************************************************************/
ChMatrixDynamic<double> H5FileInfo::GetHydrostaticStiffnessMatrix() const { //TODO cut this func
	return lin_matrix;
}

/*******************************************************************************
* H5FileInfo::GetInfAddedMassMatrix()
* returns the added mass matrix at infinite frequency
*******************************************************************************/
ChMatrixDynamic<double> H5FileInfo::GetInfAddedMassMatrix() const {
	return inf_added_mass * rho;
}

/*******************************************************************************
* H5FileInfo::GetEquilibriumCoG()
* returns cg, center of gravity of object's body
*******************************************************************************/
ChVector<> H5FileInfo::GetEquilibriumCoG() const { //TODO cut this func
	return cg;
}

/*******************************************************************************
* H5FileInfo::GetEquilibriumCoB()
* returns cb, the rotational equilibrium
*******************************************************************************/
ChVector<> H5FileInfo::GetEquilibriumCoB() const { //TODO cut this func
	return cb;
}

/*******************************************************************************
* H5FileInfo::GetRho()
* returns the density of water, rho (kg/m^3 usually)
*******************************************************************************/
double H5FileInfo::GetRho() const { //TODO cut this func
	return rho;
}

/*******************************************************************************
* H5FileInfo::GetGravity()
* returns g, gravitational acceleration, m/s^2
*******************************************************************************/
double H5FileInfo::GetGravity() const { //TODO cut this func
	return g;
}

/*******************************************************************************
* H5FileInfo::GetDisplacementVolume()
* returns displaced volume when body at equilibrium, m^3
*******************************************************************************/
double H5FileInfo::GetDisplacementVolume() const { //TODO cut this func
	return disp_vol;
}

/*******************************************************************************
* H5FileInfo::GetRIRFval()
* returns impulse response coeff for row m, column n, step s
*******************************************************************************/
double H5FileInfo::GetRIRFval(int m, int n, int s) const {
	int index = s + rirf_dims[2] * (n + m * rirf_dims[1]);
	if (index < 0 || index >= rirf_dims[0] * rirf_dims[1] * rirf_dims[2]) {
		std::cout << "out of bounds IRF\n";
		return 0;
	}
	else {
		return rirf_matrix[index] * GetRho(); // scale radiation force by rho
	}
}

/*******************************************************************************
* H5FileInfo::GetRIRFDims(int i) returns the i-th component of the dimensions of radiation_damping_matrix
* i = [0,1,2] -> [number of rows, number of columns, number of matrices]
*******************************************************************************/
int H5FileInfo::GetRIRFDims(int i) const { //TODO cut this func????
	return rirf_dims[i];
}

/*******************************************************************************
* H5FileInfo::GetRIRFTimestep() returns the RIRF's timestep (dt)
*******************************************************************************/
double H5FileInfo::GetRIRFTimestep() const { 
	return rirf_timestep;
}

//TODO: Get B(w)

/*******************************************************************************
* H5FileInfo::GetNumFreqs()
* returns number of frequencies computed
*******************************************************************************/
double H5FileInfo::GetNumFreqs() const { //TODO cut this func????
	return freq_dims[0];
}

/*******************************************************************************
* H5FileInfo::GetOmegaMin()
* returns min value of omega
*******************************************************************************/
double H5FileInfo::GetOmegaMin() const { //TODO cut this func????
	return freq_list[0];
}

/*******************************************************************************
* H5FileInfo::GetOmegaMax()
* returns max value of omega
*******************************************************************************/
double H5FileInfo::GetOmegaMax() const {
	return freq_list[freq_dims[0] - 1];
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

/*******************************************************************************
* H5FileInfo::GetRIRFdt() returns the difference in first 2 rirf_time_vector
*******************************************************************************/
double H5FileInfo::GetRIRFdt() const {
	return rirf_time_vector[1] - rirf_time_vector[0];
}

/*******************************************************************************
* H5FileInfo::GetRIRFTimeVector()
* returns the std::vector of rirf_time_vector from h5 file
*******************************************************************************/
std::vector<double> H5FileInfo::GetRIRFTimeVector() const { //TODO cut this func
	return rirf_time_vector;
}

// =============================================================================
// HydroInputs Class Definitions
// =============================================================================

/*******************************************************************************
* HydroInputs constructor
* does nothing (yet?)
*******************************************************************************/
HydroInputs::HydroInputs() {

}

// =============================================================================
// ComponentFunc Class Definitions
// =============================================================================

/*******************************************************************************
* ComponentFunc constructor
* sets pointer to ForceFunc6d member object and index for which component
* this ComponentFunc object refers to
*******************************************************************************/
ComponentFunc::ComponentFunc(ForceFunc6d* b, int i) : base(b), index(i) {	}

/*******************************************************************************
* ComponentFunc::Clone()
* required override function since ComponentFunc inherits from ChFunction
*******************************************************************************/
ComponentFunc* ComponentFunc::Clone() const  {
	return new ComponentFunc(*this);
}

/*******************************************************************************
* ComponentFunc::Get_y()
* returns the value of the index-th coordinate of the linear restoring force
* at each time
* required override function since ComponentFunc inherits from ChFunction
*******************************************************************************/
double ComponentFunc::Get_y(double x) const {
	return base->coordinateFunc(index);
}

// =============================================================================
// ForceFunc6d Class Definitions
// =============================================================================

/*******************************************************************************
* ForceFunc6d default constructor
* initializes array of ComponentFunc objects and pointers to each force/torque
*******************************************************************************/
ForceFunc6d::ForceFunc6d() : forces{ {this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5} } {
	for (unsigned i = 0; i < 6; i++) {
		force_ptrs[i] = std::shared_ptr<ComponentFunc>(forces + i, [](ComponentFunc*) {});
		// sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is
		// automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
		// them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the
		// default deletion logic to do nothing
		// Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
	}
	chrono_force = chrono_types::make_shared<ChForce>();
	chrono_torque = chrono_types::make_shared<ChForce>();
}

/*******************************************************************************
* ForceFunc6d constructor
* calls default constructor and initializes hydro force info
* from H5FileInfo
* also initializes ChBody that this force will be applied to
*******************************************************************************/
ForceFunc6d::ForceFunc6d( std::shared_ptr<ChBody> object, TestHydro* user_all_forces) : ForceFunc6d() {
	body = object;
	std::string temp = body->GetNameString(); // remove "body" from "bodyN", convert N to int, get body num
	b_num = stoi(temp.erase(0, 4));
	all_hydro_forces = user_all_forces; // TODO switch to smart pointers?
	SetForce();
	SetTorque();
	// define wave inputs here
	// TODO: switch depending on wave option (regular, regularCIC, irregular, noWaveCIC)
	//wave_amplitude = hydro_inputs.GetRegularWaveAmplitude();
	//wave_omega = hydro_inputs.GetRegularWaveOmega();
	//wave_omega_delta = file_info.GetOmegaDelta();
	//freq_index_des = (wave_omega / wave_omega_delta) - 1;

	//for (int rowEx = 0; rowEx < 6; rowEx++) {
	//	excitation_force_mag[rowEx] = file_info.GetExcitationMagInterp(rowEx, 0, freq_index_des);
	//	excitation_force_phase[rowEx] = file_info.GetExcitationPhaseInterp(rowEx, 0, freq_index_des);
	//}
}

///*******************************************************************************
//* ForceFunc6d::ComputeForceExcitationRegularFreq()
//*
//*******************************************************************************/
//ChVectorN<double, 6> ForceFunc6d::ComputeForceExcitationRegularFreq() {
//	// TODO move this check for all ComputForce type functions to just before they are called?
//	if (body->GetChTime() == previous_time_ex) {
//		return force_excitation_freq;
//	}
//	previous_time_ex = body->GetChTime();
//	for (int rowEx = 0; rowEx < 6; rowEx++) {
//		if (rowEx == 2) {
//			force_excitation_freq[rowEx] = excitation_force_mag[rowEx] * wave_amplitude * cos(wave_omega * body->GetChTime() + excitation_force_phase[rowEx]);
//		}
//		else {
//			force_excitation_freq[rowEx] = 0.0;
//		}
//	}
//	return force_excitation_freq;
//}
//
/*******************************************************************************
* ForceFunc6d::coordinateFunc
* if index is in [0,6] the corresponding vector component of the force vector
* is returned
* b_num is 1 indexed!!!!!!!!
* otherwise a warning is printed and the force is interpreted to be 0
*******************************************************************************/
double ForceFunc6d::coordinateFunc(int i) {
	return all_hydro_forces->coordinateFunc(b_num, i);
}

/*******************************************************************************
* ForceFunc6d::SetForce
* used to initialize components of force (external ChForce pointer)
*******************************************************************************/
void ForceFunc6d::SetForce() {
	chrono_force->SetF_x(force_ptrs[0]);
	chrono_force->SetF_y(force_ptrs[1]);
	chrono_force->SetF_z(force_ptrs[2]);
	body->AddForce(chrono_force);
}

/*******************************************************************************
* ForceFunc6d::SetTorque
* used to initialize components of torque (external ChForce pointer with TORQUE flag set)
*******************************************************************************/
void ForceFunc6d::SetTorque() {
	chrono_torque->SetF_x(force_ptrs[3]);
	chrono_torque->SetF_y(force_ptrs[4]);
	chrono_torque->SetF_z(force_ptrs[5]);
	chrono_torque->SetMode(ChForce::ForceType::TORQUE);
	body->AddForce(chrono_torque);
}

// =============================================================================
// TestHydro Class Definitions
// =============================================================================
TestHydro::TestHydro() {
	//std::vector<ChVectorN<double, 6>> velocity_history; //TODO vel history is gonna be tricky
	prev_time = -1;
	offset_rirf = 0;
	num_bodies = 0;
}

TestHydro::TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies, std::string h5_file_name, HydroInputs user_hydro_inputs) : TestHydro() {
	bodies = user_bodies; // 0 indexed
	num_bodies = bodies.size();
	for (int b = 0; b < num_bodies; b++) {
		file_info.emplace_back(h5_file_name, bodies[b]->GetNameString()); // set up vector of file infos for each body
	}
	hydro_inputs = user_hydro_inputs;
	// set up time vector (should be the same for each body, so just use the first always)
	rirf_time_vector = file_info[0].GetRIRFTimeVector();
	// simplify 6* num_bodies to be the system's total number of dofs, makes expressions later easier to read
	unsigned total_dofs = 6 * num_bodies;
	// resize and initialize velocity history vector to all zeros
	velocity_history.resize(file_info[0].GetRIRFDims(2) * total_dofs, 0); // resize and fill with 0s
	// resize and initialize all persistent forces to all 0s
	force_hydrostatic.resize(total_dofs, 0);
	force_radiation_damping.resize(total_dofs, 0);
	total_force.resize(total_dofs, 0);
	// set up equilibrium for entire system (each body has position and rotation equilibria 3 indicies apart)
	equilibrium.resize(total_dofs);
	for (int b = 0; b < num_bodies; b++) {
		for (int i = 0; i < 3; i++) {
			unsigned equilibrium_idx = i + 6 * b;
			equilibrium[equilibrium_idx] = file_info[b].GetEquilibriumCoG()[i];
			equilibrium[equilibrium_idx + 3ul] = file_info[b].GetEquilibriumCoB()[i]; // 3ul is just 3, used to get rid of warning on 32 vs 64 bit types
		}
	}
	for (int b = 0; b < num_bodies; b++) {
		force_per_body.emplace_back(bodies[b], this);
	}
}

// step: [0,1,...,1000] (timesteps from h5 file, one velocity per step
// c: [0,..,num_bodies-1,...,numbodies*6-1] (in order of bodies, iterates over dof for each body...3 bodies c would be [0,1,...,17])
double TestHydro::getVelHistoryAllBodies(int step, int c) const {
	//if (step * c > velocity_history.size() || c < 0 ) {
	//	// uh oh
	//	std::cout << "wrong vel_history index somewhere\nreturning 0 for it\n" << std::endl;
	//	return 0;
	//}
	//else {

	//}
	int index = c % 6;
	int b = c / 6; // 0 indexed
	return velocity_history[index + 6 * b + 6 * num_bodies * step];
}

// step: [0,1,...,1000] (0 indexed, up to the final timestep in h5 file)
// b_num: [1,2,...,total_bodies] (1 indexed!, use body number in h5 file)
// index: [0,1,2,3,4,5] (0 indexed, always 0-5 for force+torque vector indexing)
double TestHydro::setVelHistory(double val, int step, int b_num, int index) {
	//if (step * b_num * index > velocity_history.size() || (b_num > num_bodies || b_num < 1) || (index > 5 || index < 0) ) {
	//	// uh oh
	//	std::cout << "wrong vel_history index somewhere\nnot gonna set anything\n" << std::endl;
	//	return 0;
	//}
	//else {

	//}
	velocity_history[index + 6 * (b_num - 1) + 6 * num_bodies * step] = val;
	return val;
}

std::vector<double> TestHydro::ComputeForceHydrostatics() {
	std::vector<double> temp;
	unsigned total_dofs = 6 * num_bodies;
	temp.resize(total_dofs);
	for (int b = 0; b < num_bodies; b++) { // temp initialized to system current pos/rot vectors
		for (int i = 0; i < 3; i++) {
			temp[i + 6 * b] = bodies[b]->GetPos().eigen()[i];
			temp[i + 3 + 6 * b] = bodies[b]->GetRot().Q_to_Euler123()[i];
		}
	}
	// make temp be displacement vector for system
	for (int i = 0; i < total_dofs; i++) {
		temp[i] = temp[i] - equilibrium[i];
	}
	// reset force_hydrostatic to 0
	std::fill(force_hydrostatic.begin(), force_hydrostatic.end(), 0);
	// re invent matrix vector multiplication
	for (int b = 0; b < num_bodies; b++) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < total_dofs; j++) { // one of i or j needs to go to total_dofs....depends on h5 hydro matrix
				force_hydrostatic[i + 6 * b] -= ((file_info[b].GetHydrostaticStiffnessMatrix()(i,j)) * temp[j]);
			}
		}
	}
	// now handle buoyancy force....
	// one roll,pitch,yaw,buoyancy value for each body
	double* rollLeverArm = new double[num_bodies];
	double* pitchLeverArm = new double[num_bodies];
	double* yawLeverArm = new double[num_bodies];
	double* buoyancy = new double[num_bodies];
	for (int b = 0; b < num_bodies; b++) { // for each body...
		rollLeverArm[b] = force_hydrostatic[3 + b * 6]; // initialize rollLeverArm to be first rotation coord (3) for body
		pitchLeverArm[b] = force_hydrostatic[4 + b * 6]; // init pitch to be second rotation coord
		yawLeverArm[b] = force_hydrostatic[5 + b * 6]; // init yaw to be third rotation coord
		buoyancy[b] = file_info[b].GetRho() * file_info[b].GetGravity() * file_info[b].GetDisplacementVolume(); // buoyancy = rho*g*Vdisp
	}
	// TODO check buoyancy force is implemented correctly!

	for (int b = 0; b < num_bodies; b++) {
		force_hydrostatic[2 + 6 * b] += buoyancy[b];
		force_hydrostatic[3 + 6 * b] += buoyancy[b] * rollLeverArm[b];
		force_hydrostatic[4 + 6 * b] += buoyancy[b] * pitchLeverArm[b];
		force_hydrostatic[5 + 6 * b] += buoyancy[b] * yawLeverArm[b];
	}
	delete[] rollLeverArm;
	delete[] pitchLeverArm;
	delete[] yawLeverArm;
	delete[] buoyancy;
	return force_hydrostatic;
}

std::vector<double> TestHydro::ComputeForceRadiationDampingConvolutionTrapz() {
	int size = file_info[0].GetRIRFDims(2);
	// "shift" everything left 1
	offset_rirf--;
	// keep offset close to 0, avoids small chance of -overflow errors in long simulations
	if (offset_rirf < -1 * size) {
		offset_rirf += size;
	}
	int numRows = 6, numCols = 6*num_bodies;
	//int bodyNum = file_info.bodyNum;
	double* timeseries = new double[numRows * numCols * size];
	double* tmp_s = new double[numRows * size];
	// define shortcuts for accessing 1D arrays as 3D (or 2D) arrays
	// TIMESERIES is for each row in RIRF, element wise multipy velocity history by RIRF slab
#define TIMESERIES(row,col,step) timeseries[(row)*numCols*size + (col)*size + (step)]
	//TMP_S ends up being a sum over the columns of TIMESERIES (total_dofs aka LDOF
#define TMP_S(row,step) tmp_s[(row)*size + (step)]
	// set last entry as velocity
	for (int i = 0; i < 3; i++) {
		for (int b = 1; b < num_bodies + 1; b++) { // body index sucks but i think this is correct...
			setVelHistory(bodies[b-1]->GetPos_dt()[i],
				(((size + offset_rirf) % size) + size) % size, b, i);
			setVelHistory(bodies[b-1]->GetWvel_par()[i],
				(((size + offset_rirf) % size) + size) % size, b, i+3);
			//velocity_history[(((size + offset_rirf) % size) + size) % size][i] = body->GetPos_dt()[i];
			//velocity_history[(((size + offset_rirf) % size) + size) % size][i + 3] = body->GetWvel_par()[i];
		}
	}
	int vi;
	std::fill(force_radiation_damping.begin(), force_radiation_damping.end(), 0);
	//#pragma omp parallel for
	for (int row = 0; row < numRows; row++) { 
		// double fDampingCol = 0.0; // what is this?
		//#pragma omp parallel for
		for (int col = 0; col < numCols; col++) { 
			for (int st = 0; st < size; st++) {
				vi = (((st + offset_rirf) % size) + size) % size; // vi takes care of circshift function from matLab
				TIMESERIES(row, col, st) = GetRIRFval(row, col, st) * getVelHistoryAllBodies(vi,col); // col now runs thru all bodies (0->11 for 2 bodies...)
				TMP_S(row, st) = TIMESERIES(row, col, st);
				if (st > 0) {
					force_radiation_damping[row] -= (TMP_S(row, st - 1) + TMP_S(row, st)) / 2.0 * (rirf_time_vector[st] - rirf_time_vector[st - 1ull]);
				}
			}
		}
	}
	// Deallocate memory
#undef TIMESERIES
#undef TMP_S
	delete[] timeseries;
	delete[] tmp_s;
	return force_radiation_damping;
}

std::vector<double> TestHydro::ComputeForceRadiationDampingConvolutionFixed() {
	int size = file_info[0].GetRIRFDims(2);
	rirf_timestep = file_info[0].GetRIRFTimestep();
	// "shift" everything left 1
	offset_rirf--;
	// keep offset close to 0, avoids small chance of -overflow errors in long simulations
	if (offset_rirf < -1 * size) {
		offset_rirf += size;
	}
	int numRows = 6, numCols = 6 * num_bodies;
	//int bodyNum = file_info.bodyNum;
	double* timeseries = new double[numRows * numCols * size];
	double* tmp_s = new double[numRows * size];
	// define shortcuts for accessing 1D arrays as 3D (or 2D) arrays
	// TIMESERIES is for each row in RIRF, element wise multipy velocity history by RIRF slab
#define TIMESERIES(row,col,step) timeseries[(row)*numCols*size + (col)*size + (step)]
	//TMP_S ends up being a sum over the columns of TIMESERIES (total_dofs aka LDOF
#define TMP_S(row,step) tmp_s[(row)*size + (step)]
	// set last entry as velocity
	for (int i = 0; i < 3; i++) {
		for (int b = 1; b < num_bodies + 1; b++) { // body index sucks but i think this is correct...
			setVelHistory(bodies[b - 1]->GetPos_dt()[i],
				(((size + offset_rirf) % size) + size) % size, b, i);
			setVelHistory(bodies[b - 1]->GetWvel_par()[i],
				(((size + offset_rirf) % size) + size) % size, b, i + 3);
			//velocity_history[(((size + offset_rirf) % size) + size) % size][i] = body->GetPos_dt()[i];
			//velocity_history[(((size + offset_rirf) % size) + size) % size][i + 3] = body->GetWvel_par()[i];
		}
	}
	int vi;
	std::fill(force_radiation_damping.begin(), force_radiation_damping.end(), 0);
	//#pragma omp parallel for
	for (int row = 0; row < numRows; row++) {
		// double fDampingCol = 0.0; // what is this?
		//#pragma omp parallel for
		sumVelHistoryAndRIRF = 0.0;
		for (int col = 0; col < numCols; col++) {
			for (int st = 0; st < size; st++) {
				vi = (((st + offset_rirf) % size) + size) % size; // vi takes care of circshift function from matLab
				TIMESERIES(row, col, st) = GetRIRFval(row, col, st) * getVelHistoryAllBodies(vi, col); // col now runs thru all bodies (0->11 for 2 bodies...)
				TMP_S(row, st) = TIMESERIES(row, col, st); //TODO: rename TIMESERIES and TMP_S
				sumVelHistoryAndRIRF += TMP_S(row, st);
			}
		}
		force_radiation_damping[row] -= sumVelHistoryAndRIRF * rirf_timestep;
	}
	// Deallocate memory
#undef TIMESERIES
#undef TMP_S
	delete[] timeseries;
	delete[] tmp_s;
	return force_radiation_damping;
}

double TestHydro::GetRIRFval(int row, int col, int st) {
	int b = col / 6; // 0 indexed
	int c = col % 6;
	return file_info[b].GetRIRFval(row, c, st);
}

// compute all the forces here, b is 1 indexed
double TestHydro::coordinateFunc(int b, int i) { // b_num from ForceFunc6d is 1 indexed
	unsigned body_num_offset = 6 * (b - 1);
	if (i < 0 || i > 5) {
		std::cout << "wrong index somewhere\nsetting coordinateFunc to 0" << std::endl;
		return 0;
	}
	// check prev_time here and only here
	// if forces have been computed for this time already, return the computed total force
	if (bodies[0]->GetChTime() == prev_time) {
		return total_force[body_num_offset + i];
	}
	// update current time and total_force for this step
	prev_time = bodies[0]->GetChTime();
	// call all compute force functions
	ComputeForceHydrostatics();
	//ComputeForceRadiationDampingConvolutionTrapz();
	ComputeForceRadiationDampingConvolutionFixed();

	// sum all forces element by element
	unsigned total_dofs = 6 * num_bodies;
	for (int j = 0; j < total_dofs; j++) {
		total_force[j] = force_hydrostatic[j] + force_radiation_damping[j];
	}
	return total_force[body_num_offset + i];
}

// =============================================================================
// ChLoadAddedMass Class Definitions
// =============================================================================

/*******************************************************************************
* constructorHelper
* for use in ChLoadAddedMass constructor, to conver a vector of shared_ptrs to 
* ChBody type objects to a vector of shared_ptrs to ChLoadable type objects
* in order to use ChLoadCustomMultiple's constructor for a vector of ChLoadable
* shared_ptrs
*******************************************************************************/
//std::vector<std::shared_ptr<ChLoadable>> constructorHelper(std::vector<std::shared_ptr<ChBody>>& bodies) {
//	std::vector<std::shared_ptr<ChLoadable>> re(bodies.size());
//
//	std::transform(bodies.begin(), bodies.end(), re.begin(), [](const std::shared_ptr<ChBody>& p) { return std::static_pointer_cast<ChLoadable>(p); });
//
//	return re;
//}

/*******************************************************************************
* ChLoadAddedMass constructor
* initializes body to have load applied to and added mass matrix from h5 file object
*******************************************************************************/
//ChLoadAddedMass::ChLoadAddedMass(const H5FileInfo& file, 
//	std::vector<std::shared_ptr<ChBody>>& bodies) 
//	: ChLoadCustomMultiple(constructorHelper(bodies)) { ///< calls ChLoadCustomMultiple to link loads to bodies
//	inf_added_mass_J = file.GetInfAddedMassMatrix(); //TODO switch all uses of H5FileInfo object to be like this, instead of copying the object each time?
//
//	std::ofstream myfile;
//	//myfile.open("C:\\code\\chrono_hydro_dev\\HydroChrono_build\\Release\\inf_added_mass_J.txt");
//	myfile.open("inf_added_mass_J_txt");
//	myfile << inf_added_mass_J << "\n";
//	myfile.close();
//}
/*******************************************************************************
* ChLoadAddedMass constructor
* initializes body to have load applied to and added mass matrix from h5 file object
*******************************************************************************/
//ChLoadAddedMass::ChLoadAddedMass(const ChMatrixDynamic<>& addedMassMatrix,
//	std::vector<std::shared_ptr<ChBody>>& bodies)
//	: ChLoadCustomMultiple(constructorHelper(bodies)) { ///< calls ChLoadCustomMultiple to link loads to bodies
//	inf_added_mass_J = addedMassMatrix; //TODO switch all uses of H5FileInfo object to be like this, instead of copying the object each time?
//
//	std::ofstream myfile;
//	//myfile.open("C:\\code\\chrono_hydro_dev\\HydroChrono_build\\Release\\inf_added_mass_J.txt");
//	myfile.open("inf_added_mass_J_txt");
//	myfile << inf_added_mass_J << "\n";
//	myfile.close();
//}
/*******************************************************************************
* ChLoadAddedMass::ComputeJacobian()
* Computes Jacobian for load, in this case just the mass matrix is initialized
* as the added mass matrix
*******************************************************************************/
//void ChLoadAddedMass::ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
//	ChStateDelta* state_w,  ///< state speed to evaluate jacobians
//	ChMatrixRef mK,         ///< result dQ/dx
//	ChMatrixRef mR,         ///< result dQ/dv
//	ChMatrixRef mM          ///< result dQ/da
//) {
//	//set mass matrix here
//
//	jacobians->M = inf_added_mass_J;
//
//	ChMatrixDynamic<double> massmat = jacobians->M;
//
//	std::ofstream myfile2;
//	//myfile2.open("C:\\code\\chrono_hydro_dev\\HydroChrono_build\\Release\\massmat.txt");
//	myfile2.open("massmat.txt");
//	myfile2 << massmat << "\n";
//	myfile2.close();
//
//	// R gyroscopic damping matrix terms (6x6)
//	// 0 for added mass
//	jacobians->R.setZero();
//
//	// K inertial stiffness matrix terms (6x6)
//	// 0 for added mass
//	jacobians->K.setZero();
//}

/*******************************************************************************
* ChLoadAddedMass::LoadIntLoadResidual_Mv()
* Computes LoadIntLoadResidual_Mv for vector w, const c, and vector R
* Note R here is vector, and is not R gyroscopic damping matrix from ComputeJacobian
*******************************************************************************/
//void ChLoadAddedMass::LoadIntLoadResidual_Mv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) {
//	if (!this->jacobians)
//		return;
//
//	//if (!loadable->IsSubBlockActive(0))
//	//	return;
//
//	// R+=c*M*a
//	// segment gives the chunk of vector starting at the first argument, and going for as many elements as the second argument...
//	// in this case, segment gets the 3vector starting at the 0th DOF's offset (ie 0)
//	//R.segment(loadable->GetSubBlockOffset(0), 3) += c * (this->mass * (a_x + chrono::Vcross(a_w, this->c_m))).eigen();
//	// in this case, segment gets the 3vector starting at the 0th DOF's + 3 offset (ie 3)
//	//R.segment(loadable->GetSubBlockOffset(0) + 3, 3) += c * (this->mass * chrono::Vcross(this->c_m, a_x) + this->I * a_w).eigen();
//	// since R is a vector, we can probably just do R += C*M*a with no need to separate w into a_x and a_w above
//	R += c * jacobians->M * w;
//}

// =============================================================================
// LoadAllHydroForces Class Definitions
// =============================================================================
/*******************************************************************************
*
*******************************************************************************/
//LoadAllHydroForces::LoadAllHydroForces(std::vector<std::shared_ptr<ChBody>> objects, std::string file/*, std::string bodyName*/, HydroInputs user_hydro_inputs)/* :
//	sys_file_info(file, bodyName), hydro_force(sys_file_info, object, user_hydro_inputs)*/ {
//	
//	//TestHydro(std::vector<ForceFunc6d>&fpb, HydroInputs users_hydro_inputs);
//
//	for (int i = 0; i < objects.size(); i++) { 
//		
//	}
//	//hydro_force.SetForce();
//	//hydro_force.SetTorque(); // called in constructor for 6d class now
//}
