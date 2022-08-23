#include "hydro_forces.h"

#include <algorithm>
#include <numeric> // std::accumulate

// =============================================================================
// H5FileInfo Class Definitions
// =============================================================================

H5FileInfo::H5FileInfo() {}

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

	//if (printed == false) {
	//	std::ofstream H5FileOut;
	//	H5FileOut.open("C:\\code\\HydroChrono_build\\Release\\results\\rm3\\debugging\\H5FileOut.txt");
	//	H5FileOut << Name << "\n";
	//	H5FileOut << bodyName << "\n";
	//	printed = true;
	//	H5FileOut.close();
	//}

}

/*******************************************************************************
* H5FileInfo::operator = (H5FileInfo& rhs) 
* defines = operator, sets the left object = to right object (*this = rhs)
* returns *this
*******************************************************************************/
H5FileInfo& H5FileInfo::operator = (H5FileInfo& rhs) {
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
	mag_dims = rhs.mag_dims;
	excitation_re_matrix = rhs.excitation_re_matrix;
	re_dims = rhs.re_dims;
	excitation_im_matrix = rhs.excitation_im_matrix;
	im_dims = rhs.im_dims;
	h5_file_name = rhs.h5_file_name;
	bodyName = rhs.bodyName;
	return *this;
}

/*******************************************************************************
* H5FileInfo copy constructor (H5FileInfo& rhs)
* defines basic copy constructor using the = operator
*******************************************************************************/
H5FileInfo::H5FileInfo(H5FileInfo& old) {
	*this = old;
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

	Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/mag", excitation_mag_matrix, mag_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/re", excitation_re_matrix, re_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/excitation/im", excitation_im_matrix, im_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/radiation_damping/impulse_response_fun/K", rirf_matrix, rirf_dims);
	Init3D(userH5File, bodyName + "/hydro_coeffs/radiation_damping/all", radiation_damping_matrix, Bw_dims);
	
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
//double H5FileInfo::GetNumFreqs() const { //TODO cut this func????
//	return freq_dims[0];
//}

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
//double H5FileInfo::GetOmegaMax() const {
//	return freq_list[freq_dims[0] - 1];
//}

/*******************************************************************************
* H5FileInfo::GetOmegaDelta()
* returns omega step size
*******************************************************************************/
//double H5FileInfo::GetOmegaDelta() const {
//	return GetOmegaMax() / GetNumFreqs();
//}

/*******************************************************************************
* H5FileInfo::GetExcitationMagValue()
* returns excitation magnitudes for row i, column j, frequency ix k
*******************************************************************************/
//double H5FileInfo::GetExcitationMagValue(int i, int j, int k) const {
//	int indexExMag = k + excitation_mag_dims[2] * i;
//	return excitation_mag_matrix[indexExMag] * rho * g;
//}

/*******************************************************************************
* H5FileInfo::GetExcitationMagInterp()
* returns excitation magnitudes for row i, column j, frequency ix k
*******************************************************************************/
//double H5FileInfo::GetExcitationMagInterp(int i, int j, double freq_index_des) const {
//	double freq_interp_val = freq_index_des - floor(freq_index_des);
//	double excitationMagFloor = GetExcitationMagValue(i, j, floor(freq_index_des));
//	double excitationMagCeil = GetExcitationMagValue(i, j, floor(freq_index_des) +1);
//	double excitationMag = (freq_interp_val * (excitationMagCeil - excitationMagFloor)) + excitationMagFloor;
//
//	return excitationMag;
//}

/*******************************************************************************
* H5FileInfo::GetExcitationPhaseValue()
* returns excitation phases for row i, column j, frequency k
*******************************************************************************/
//double H5FileInfo::GetExcitationPhaseValue(int i, int j, int k) const {
//	int indexExPhase = k + excitation_phase_dims[2] * i;
//	return excitation_phase_matrix[indexExPhase];
//}

/*******************************************************************************
* H5FileInfo::GetExcitationPhaseInterp()
* returns excitation phases for row i, column j, frequency ix k
*******************************************************************************/
//double H5FileInfo::GetExcitationPhaseInterp(int i, int j, double freq_index_des) const {
//	double freq_interp_val = freq_index_des - floor(freq_index_des);
//	double excitationPhaseFloor = GetExcitationPhaseValue(i, j, floor(freq_index_des));
//	double excitationPhaseCeil = GetExcitationPhaseValue(i, j, floor(freq_index_des) + 1);
//	double excitationPhase = (freq_interp_val * (excitationPhaseCeil - excitationPhaseFloor)) + excitationPhaseFloor;
//
//	return excitationPhase;
//}

/*******************************************************************************
* H5FileInfo::GetRIRFdt() returns the difference in first 2 rirf_time_vector
*******************************************************************************/
//double H5FileInfo::GetRIRFdt() const {
//	return rirf_time_vector[1] - rirf_time_vector[0];
//}

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
* ComponentFunc default constructor, dont use
* sets pointer to ForceFunc6d member object to null and index to 6 (invalid)
* this ComponentFunc object refers to
*******************************************************************************/
ComponentFunc::ComponentFunc() {
	base = NULL;
	index = 6;
}

/*******************************************************************************
* ComponentFunc constructor
* sets pointer to ForceFunc6d member object and index for which component
* this ComponentFunc object refers to
*******************************************************************************/
ComponentFunc::ComponentFunc(ForceFunc6d* b, int i) : base(b), index(i) { }

/*******************************************************************************
* ComponentFunc::Clone()
* required override function since ComponentFunc inherits from ChFunction
*******************************************************************************/
ComponentFunc* ComponentFunc::Clone() const  {
	return new ComponentFunc(*this);
}

/*******************************************************************************
* ComponentFunc copy constructor
* copy constructor works like defualt
*******************************************************************************/
ComponentFunc::ComponentFunc(const ComponentFunc& old) {
	base = old.base;
	index = old.index;
}

/*******************************************************************************
* ComponentFunc::Get_y()
* returns the value of the index-th coordinate of the linear restoring force
* at each time
* required override function since ComponentFunc inherits from ChFunction
*******************************************************************************/
double ComponentFunc::Get_y(double x) const {
	if (base == NULL) {
		std::cout << "base == Null!" << std::endl;
		return 0;
	}
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
	b_num = stoi(temp.erase(0, 4)); // 1 indexed
	all_hydro_forces = user_all_forces; // TODO switch to smart pointers? does this use = ?
	if (all_hydro_forces == NULL) {
		std::cout << "all hydro forces null " << std::endl;
	}
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

/*******************************************************************************
* ForceFunc6d copy constructor
* calls
*******************************************************************************/
ForceFunc6d::ForceFunc6d(const ForceFunc6d& old) : 
forces{ {this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5} } {
	for (unsigned i = 0; i < 6; i++) {
		force_ptrs[i] = std::shared_ptr<ComponentFunc>(forces + i, [](ComponentFunc*) {});
		// sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is
		// automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
		// them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the
		// default deletion logic to do nothing
		// Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
	}
	chrono_force = old.chrono_force;
	chrono_torque = old.chrono_torque;
	body = old.body;
	b_num = old.b_num;
	all_hydro_forces = old.all_hydro_forces;
	SetForce();
	SetTorque();
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
	// b_num is 1 indexed?
	if (i >= 6 || i < 0) {
		std::cout << "wrong index force func 6d" << std::endl;
		return 0;
	}
	return all_hydro_forces->coordinateFunc(b_num, i);
}

/*******************************************************************************
* ForceFunc6d::SetForce
* used to initialize components of force (external ChForce pointer)
*******************************************************************************/
void ForceFunc6d::SetForce() {
	if (chrono_force == NULL || body == NULL) {
		std::cout << "set force null issue" << std::endl;
	}
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
	if (chrono_torque == NULL || body == NULL) {
		std::cout << "set torque null issue" << std::endl;
	}
	chrono_torque->SetF_x(force_ptrs[3]);
	chrono_torque->SetF_y(force_ptrs[4]);
	chrono_torque->SetF_z(force_ptrs[5]);
	chrono_torque->SetMode(ChForce::ForceType::TORQUE);
	body->AddForce(chrono_torque);
}

// =============================================================================
// TestHydro Class Definitions
// =============================================================================
/*******************************************************************************
* TestHydro::TestHydro()
* default constructor to initialize just a few variables
* called from actual constructor
*******************************************************************************/
TestHydro::TestHydro() {
	//std::vector<ChVectorN<double, 6>> velocity_history; //TODO vel history is gonna be tricky
	prev_time = -1;
	offset_rirf = 0;
	num_bodies = 0;
}

/*******************************************************************************
* TestHydro::TestHydro(user_bodies, h5_file_name, user_hydro_inputs)
* main constructor for TestHydro class, sets up vector of bodies, h5 file info,
* and hydro inputs
* also initializes many persistent variables for force calculations
* calls default constructor 
*******************************************************************************/
TestHydro::TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies, std::string h5_file_name, HydroInputs user_hydro_inputs) : TestHydro() {
	bodies = user_bodies; // 0 indexed
	num_bodies = bodies.size();
	for (int b = 0; b < num_bodies; b++) {
		file_info.emplace_back(h5_file_name, bodies[b]->GetNameString()); // set up vector of file infos for each body
	}
	//hydro_inputs = user_hydro_inputs;
	// set up time vector (should be the same for each body, so just use the first always)
	rirf_time_vector = file_info[0].GetRIRFTimeVector();
	rirf_timestep = rirf_time_vector[1] - rirf_time_vector[0]; // TODO is this the same for all bodies?
	// simplify 6* num_bodies to be the system's total number of dofs, makes expressions later easier to read
	unsigned total_dofs = 6 * num_bodies;
	// resize and initialize velocity history vector to all zeros
	velocity_history.resize(file_info[0].GetRIRFDims(2) * total_dofs, 0); // resize and fill with 0s
	// resize and initialize all persistent forces to all 0s
	force_hydrostatic.resize(total_dofs, 0.0);
	force_radiation_damping.resize(total_dofs, 0.0);
	total_force.resize(total_dofs, 0.0);
	// set up equilibrium for entire system (each body has position and rotation equilibria 3 indicies apart)
	equilibrium.resize(total_dofs, 0);
	cb_minus_cg.resize(3 * num_bodies, 0);
	for (int b = 0; b < num_bodies; b++) {
		for (int i = 0; i < 3; i++) {
			unsigned equilibrium_idx = i + 6 * b;
			unsigned c_idx = i + 3 * b;
			equilibrium[equilibrium_idx] = file_info[b].cg[i]; // positional equilib is cg, leave rotational bit 0
			cb_minus_cg[c_idx] = file_info[b].cb[i] - file_info[b].cg[i];
		}
	}

	for (int b = 0; b < num_bodies; b++) {
		if (this == NULL) {
			std::cout << "woops" << std::endl;
		}
		force_per_body.emplace_back(bodies[b], this);
	}

	my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();
	my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(file_info, bodies);
	bodies[0]->GetSystem()->Add(my_loadcontainer);
	my_loadcontainer->Add(my_loadbodyinertia);
}

/*******************************************************************************
* TestHydro::getVelHistoryAllBodies(int step, int c) const
* finds and returns the component of velocity history for given step and c (column)
* step: [0,1,...,1000] (timesteps from h5 file, one velocity per step
* c: [0,..,num_bodies-1,...,numbodies*6-1] (in order of bodies, iterates over 
*    dof for each body...3 bodies c would be [0,1,...,17])
*******************************************************************************/
double TestHydro::getVelHistoryAllBodies(int step, int c) const {
	if (step < 0 || step >= file_info[0].GetRIRFDims(2) || c < 0 || c >= num_bodies * 6) {
		std::cout << "wrong vel history index " << std::endl;
		return 0;
	}
	int index = c % 6;
	int b = c / 6; // 0 indexed
	if (index + 6 * b + 6 * num_bodies * step >= num_bodies * 6 * file_info[0].GetRIRFDims(2) || index + 6 * b + 6 * num_bodies * step < 0) {
		std::cout << "bad vel history math" << std::endl;
		return 0; 
	}
	return velocity_history[index + 6 * b + 6 * num_bodies * step];
}

/*******************************************************************************
* TestHydro::setVelHistoryAllBodies(double val, int step, int b_num, int index) 
* sets velocity history for step, b_num (body number) and index (dof) to the given val
* val: value to set the requested element to
* step: [0,1,...,1000] (0 indexed, up to the final timestep in h5 file)
* b_num: [1,2,...,total_bodies] (1 indexed!, use body number in h5 file)
* index: [0,1,2,3,4,5] (0 indexed, always 0-5 for force+torque vector indexing)
*******************************************************************************/
double TestHydro::setVelHistory(double val, int step, int b_num, int index) {
	if (step < 0 || step >= file_info[0].GetRIRFDims(2) || b_num < 1 || b_num > num_bodies || index < 0 || index >= 6) {
		std::cout << "bad set vel history indexing" << std::endl;
		return 0;
	}
	if (index + 6 * (b_num - 1) + 6 * num_bodies * step < 0 || index + 6 * (b_num - 1) + 6 * num_bodies * step >= num_bodies * 6 * file_info[0].GetRIRFDims(2)) {
		std::cout << "bad set vel history math" << std::endl;
		return 0;
	}
	velocity_history[index + 6 * (b_num - 1) + 6 * num_bodies * step] = val;
	return val;
}

/*******************************************************************************
* TestHydro::ComputeForceHydrostatics() 
* computes the 6N dimensional Hydrostatic stiffness force
*******************************************************************************/
std::vector<double> TestHydro::ComputeForceHydrostatics() {
	std::vector<double> position, displacement;
	unsigned total_dofs = 6 * num_bodies;
	position.resize(total_dofs,0);
	displacement.resize(total_dofs, 0);

	std::ofstream GetPosOut;
	GetPosOut.open("results/rm3/debugging/GetPosOut.txt");
	// orientation initialized to system current pos/rot vectors
	for (int b = 0; b < num_bodies; b++) { 
		for (int i = 0; i < 3; i++) {
			unsigned b_offset = 6 * b;
			if (b_offset + i + 3 > total_dofs || b_offset + i < 0) {
				std::cout << "temp index in hydrostatic force is bad " << std::endl;
			}

			GetPosOut << b << "\t" << i << "\t" << bodies[b]->GetPos().eigen()[i] << "\t\t" << equilibrium[i + (6 * b)] << "\t\t" << bodies[b]->GetPos().eigen()[i] - equilibrium[i + (6 * b)] << "\n";
			
			//double linearPosition = bodies[b]->GetPos().eigen()[i];
			//double rotationalPosition = bodies[b]->GetRot().Q_to_Euler123()[i];
			position[i + b_offset] = bodies[b]->GetPos().eigen()[i];
			//position[i + 3 + b_offset] = rotationalPosition;
		}
	}
	GetPosOut.close();

	//std::ofstream positionOut;
	//positionOut.open("results/rm3/debugging/position.txt");
	//for (int b = 0; b < num_bodies; b++) { // for each body...
	//	for (int i = 0; i < 6; i++) {
	//		positionOut << position[i + (6 * b)] << "\n";
	//	}
	//}
	//positionOut.close();

	//std::ofstream equilibriumOut;
	//equilibriumOut.open("results/rm3/debugging/equilibrium.txt");
	//for (int b = 0; b < num_bodies; b++) { // for each body...
	//	for (int i = 0; i < 6; i++) {
	//		equilibriumOut << equilibrium[i + (6 * b)] << "\n";
	//	}
	//}
	//equilibriumOut.close();

	// make displacement vector for system
	for (int i = 0; i < total_dofs; i++) {
		displacement[i] = position[i] - equilibrium[i];
	}

	// reset force_hydrostatic to 0
	std::fill(force_hydrostatic.begin(), force_hydrostatic.end(), 0.0);

	std::ofstream dispOut;
	dispOut.open("results/rm3/debugging/dispOut.txt");
	for (int i = 0; i < total_dofs; i++) {
		dispOut << displacement[i] << "\n";
	}

	std::ofstream fHSOut0;
	fHSOut0.open("results/rm3/debugging/fHSOut0.txt");
	for (int i = 0; i < total_dofs; i++) {
		fHSOut0 << force_hydrostatic[i] << "\n";
	}

	std::ofstream KHS0Out;
	KHS0Out.open("results/rm3/debugging/KHS0Out.txt");
	KHS0Out << file_info[0].lin_matrix << "\n";

	std::ofstream KHS1Out;
	KHS1Out.open("results/rm3/debugging/KHS1Out.txt");
	KHS1Out << file_info[1].lin_matrix << "\n";
	// re invent matrix vector multiplication
	for (int b = 0; b < num_bodies; b++) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
	//			//KHSOut << i << " " << j << " " << file_info[b].GetHydrostaticStiffness(i, j) << "\n";
				force_hydrostatic[i + (6 * b)] -= ((file_info[b].GetHydrostaticStiffness(i, j)) * displacement[i + (6 * b)]);
			}
		}
	}

	std::ofstream fHSOut1;
	fHSOut1.open("results/rm3/debugging/fHSOut1.txt");
	for (int b = 0; b < num_bodies; b++) {
		for (int i = 0; i < 6; i++) {
			fHSOut1 << force_hydrostatic[i + (6*b)] << "\t\t" << (1000 * 9.81 * file_info[b].disp_vol) << "\t\t" << -(bodies[b]->GetMass() * 9.81) << "\n";
		}
	}


	////// now handle buoyancy force....
	//std::ofstream displacementOut;
	//displacementOut.open("results/rm3/debugging/displacement.txt");
	//for (int b = 0; b < num_bodies; b++) { // for each body...
	//	for (int i = 0; i < 6; i++) {
	//		displacementOut << displacement[i + (6 * b)] << "\n";
	//	}
	//}
	//displacementOut.close();

	////// now handle buoyancy force....
	//std::ofstream force_hydrostaticOut;
	//force_hydrostaticOut.open("results/rm3/debugging/force_hydrostatic.txt");
	//for (int b = 0; b < num_bodies; b++) { // for each body...
	//	for (int i = 0; i < 6; i++) {
	//		force_hydrostaticOut << force_hydrostatic[i + (6 * b)] << "\n";
	//	}
	//}
	//force_hydrostaticOut.close();

	
	
	//std::ofstream buoyancyOut;
	//buoyancyOut.open("results/rm3/debugging/buoyancy.txt");
	//for (int b = 0; b < num_bodies; b++) { // for each body...
	//	buoyancy[b] = file_info[b].rho * file_info[b].g * file_info[b].disp_vol; // buoyancy = rho*g*Vdisp
	//	buoyancyOut << file_info[b].rho << "\n";
	//	buoyancyOut << file_info[b].g << "\n";
	//	buoyancyOut << file_info[b].disp_vol << "\n";
	//	buoyancyOut << buoyancy[b] << "\n\n";
	//}
	//buoyancyOut.close();
	// 
	
	// now handle buoyancy force....
	assert(num_bodies > 0);
	double* buoyancy = new double[num_bodies];
	double* weight = new double[num_bodies];
	// add vertical buoyancy for each body, and add (0,0,buoyancy)x(cb-cg) to torque for each body (simplified)
	for (int b = 0; b < num_bodies; b++) {
		buoyancy[b] = 1000.0 * 9.81 * file_info[b].disp_vol; // buoyancy = rho*g*Vdisp
		weight[b] = (bodies[b]->GetMass() * 9.81);
		unsigned b_offset = 6 * b;
		force_hydrostatic[2 + b_offset] += buoyancy[b];// -weight[b];////buoyancy[b];// ; // add regular z direction buoyancy force
		force_hydrostatic[2 + b_offset] -= weight[b];
		//force_hydrostatic[3 + b_offset] += -1 * buoyancy[b] * cb_minus_cg[1]; // roll part of cross product simplified
		//force_hydrostatic[4 + b_offset] += buoyancy[b] * cb_minus_cg[0]; // pitch part of cross product simplified
	}

	std::ofstream fHSOut2;
	fHSOut2.open("results/rm3/debugging/fHSOut2.txt");
	for (int i = 0; i < total_dofs; i++) {
		fHSOut2 << force_hydrostatic[i] << "\n";
	}

	delete[] buoyancy;
	delete[] weight;
	return force_hydrostatic;
}

/*******************************************************************************
* TestHydro::ComputeForceRadiationDampingConv()
* computes the 6N dimensional Radiation Damping force with convolution history
*******************************************************************************/
std::vector<double> TestHydro::ComputeForceRadiationDampingConv() {
	int size = file_info[0].GetRIRFDims(2);
	int nDoF = 6;
	// "shift" everything left 1
	offset_rirf--;
	// keep offset close to 0, avoids small chance of -overflow errors in long simulations
	if (offset_rirf < -1 * size) {
		offset_rirf += size;
	}
	int numRows = nDoF * num_bodies;
	int numCols = nDoF * num_bodies;
	assert(numRows * size > 0 && numCols > 0);
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
		}
	}
	int vi;
	std::fill(force_radiation_damping.begin(), force_radiation_damping.end(), 0);
	//#pragma omp parallel for
	if (convTrapz == true){
	// convolution integral using trapezoidal rule
		for (int row = 0; row < numRows; row++) { // row goes to 6N
			for (int st = 0; st < size; st++) {
				vi = (((st + offset_rirf) % size) + size) % size; // vi takes care of circshift function from matLab
				TMP_S(row, st) = 0;
				for (int col = 0; col < numCols; col++) { // numCols goes to 6N
					// multiply rirf by velocity history for each step and row (0,...,6N), store product in TIMESERIES
					TIMESERIES(row, col, st) = GetRIRFval(row, col, st) * getVelHistoryAllBodies(vi, col); 
					// TMP_S is the sum over col (sum the effects of all radiating dofs (LDOF) for each time and motion dof)
					TMP_S(row, st) += TIMESERIES(row, col, st);
				}
				if (st > 0) {
					// integrate TMP_S
					force_radiation_damping[row] += (TMP_S(row, st - 1) + TMP_S(row, st)) / 2.0 * (rirf_time_vector[st] - rirf_time_vector[st - 1]);
				}
			}
		}
	}
	else { // TODO fix this for force_radiation_damping to go over col not row!
		// convolution integral assuming fixed dt
		for (int row = 0; row < numRows; row++) {
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
	}
	// Deallocate memory
#undef TIMESERIES
#undef TMP_S
	delete[] timeseries;
	delete[] tmp_s;

	std::ofstream fRadOut("results/rm3/debugging/fRadOut.txt");
	for (int i = 0; i < numCols; i++) {
		fRadOut << force_radiation_damping[i] << std::endl;
	}
	fRadOut.close();

	return force_radiation_damping;
}

//std::vector<double> TestHydro::ComputeForceExcitation() {
//	// do math
//	// return vector for that force
//	force_excitation.resize(6, 0);
//	return force_excitation;
//}

/*******************************************************************************
* TestHydro::GetRIRFval(int row, int col, int st)
* returns the rirf value from the correct body given the row, step, and index 
* row: encodes the body number and dof index [0,...,5,...6N-1] for rows of RIRF
* col: col in RIRF matrix [0,...,5,...6N-1]
* st: which step in rirf ranges usually [0,...1000]
* c: gets just the index aka dof from col
* b: gets body num from col
*******************************************************************************/
double TestHydro::GetRIRFval(int row, int col, int st) {
	if (row < 0 || row >= 6 * num_bodies || col < 0 || col >= 6 * num_bodies || st < 0 || st >= file_info[0].GetRIRFDims(2)) {
		std::cout << "rirfval index bad from testhydro" << std::endl;
		return 0;
	}
	int b = row / 6; // 0 indexed, which body to get matrix info from
	int c = col % 6; // which dof across column, 0,..,11 for 2 bodies, 0,...,6N-1 for N
	int r = row % 6; // which dof 0,..,5 in individual body RIRF matrix
	return file_info[b].GetRIRFval(r, col, st);
}

/*******************************************************************************
* TestHydro::GetRIRFval(int row, int col, int st)
* ****IMPORTANT b is 1 indexed here
* computes all forces or returns saved force if it's already been calculated this
* timestep
* b: body number, it is 1 indexed here since it comes from ForcFunc6d
* i: index or Degree of Freedom (dof) ranges [0,...5]
* calls computeForce type functions
*******************************************************************************/
double TestHydro::coordinateFunc(int b, int i) { // b_num from ForceFunc6d is 1 indexed, TODO: make all b_num 0 indexed
	unsigned body_num_offset = 6 * (b - 1);
	if (i < 0 || i > 5 || b < 1 || b > num_bodies) {
		std::cout << "wrong index somewhere\nsetting coordinateFunc to 0" << std::endl;
		return 0;
	}
	// check prev_time here and only here
	// if forces have been computed for this time already, return the computed total force
	if (bodies[0] == NULL) {
		std::cout << "bodies empty" << std::endl;
		return 0;
	}
	if (bodies[0]->GetChTime() == prev_time) {
		return total_force[body_num_offset + i];
	}
	// update current time and total_force for this step
	prev_time = bodies[0]->GetChTime();
	// call all compute force functions
	ComputeForceHydrostatics();
	convTrapz = true; // use trapeziodal rule or assume fixed dt.
	ComputeForceRadiationDampingConv();

	unsigned total_dofs = 6 * num_bodies;

	// sum all forces element by element
	for (int j = 0; j < total_dofs; j++) {
		total_force[j] = force_hydrostatic[j] - force_radiation_damping[j];
	}
	if (body_num_offset + i < 0 || body_num_offset >= total_dofs) {
		std::cout << "total force accessing out of bounds" << std::endl;
	}

	//if (printed == false) {
	//	std::ofstream hsdebug;
	//	hsdebug.open("results/rm3/debugging/force_hydrostatic.txt");
	//	for (int j = 0; j < total_dofs; j++) {
	//		hsdebug << total_force[j] << "\n";
	//	}
	//	hsdebug << "\n";
	//	printed = true;
	//	hsdebug.close();
	//}

	return total_force[body_num_offset + i];
}

// =============================================================================
// ChLoadAddedMass Class Definitions
// =============================================================================

/*******************************************************************************
* constructorHelper
* for use in ChLoadAddedMass constructor, to convert a vector of shared_ptrs to 
* ChBody type objects to a vector of shared_ptrs to ChLoadable type objects
* in order to use ChLoadCustomMultiple's constructor for a vector of ChLoadable
* shared_ptrs
* TLDR: this function manually casts std::vector<std::shared_ptr<ChBody>> to a 
* std::vector<std::shared_ptr<ChLoadable>> so the ChLoadCustomMultiple constructor 
* can attach the load to the body automatically...since ChBody is a ChLoadable
* object, this cast is usually implicitly defined-the reason it isn't here is 
* because of the vector and shared_ptr parts making it a bit too complicated
* 
* re is just the object to be returned
*******************************************************************************/
std::vector<std::shared_ptr<ChLoadable>> constructorHelper(std::vector<std::shared_ptr<ChBody>>& bodies) {
	std::vector<std::shared_ptr<ChLoadable>> re(bodies.size());

	std::transform(bodies.begin(), bodies.end(), re.begin(), [](const std::shared_ptr<ChBody>& p) { return std::static_pointer_cast<ChLoadable>(p); });

	return re;
}

/*******************************************************************************
* ChLoadAddedMass::AssembleSystemAddedMassMat()
* initializes infinite_added_mass matrix by 
*******************************************************************************/
void ChLoadAddedMass::AssembleSystemAddedMassMat() {
	infinite_added_mass.setZero(6 * nBodies, 6 * nBodies);
	for (int i = 0; i < nBodies; i++) {
		infinite_added_mass.block(i * 6, 0, 6, nBodies * 6) = h5_body_data[i].GetInfAddedMassMatrix();
	}
}

/*******************************************************************************
* ChLoadAddedMass constructor
* initializes body to have load applied to and added mass matrix from h5 file object
*******************************************************************************/
ChLoadAddedMass::ChLoadAddedMass(const std::vector<H5FileInfo>& user_h5_body_data,
								 std::vector<std::shared_ptr<ChBody>>& bodies)
								     : ChLoadCustomMultiple(constructorHelper(bodies)) { ///< calls ChLoadCustomMultiple to link loads to bodies
	nBodies = bodies.size();
	h5_body_data = user_h5_body_data;
	AssembleSystemAddedMassMat();

	//ChMatrixDynamic<double> massmat = infinite_added_mass;
	//std::ofstream myfile2;
	//myfile2.open("results/rm3/debugging/massmat1.txt");
	//myfile2 << massmat << std::endl;
	//myfile2.close();
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
	jacobians->M = infinite_added_mass;

	//ChMatrixDynamic<double> massmat = jacobians->M;
	//std::ofstream myfile2;
	//myfile2.open("results/rm3/debugging/massmat1.txt");
	//myfile2 << massmat << std::endl;
	//myfile2.close();

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

	//if (!loadable->IsSubBlockActive(0))
	//	return;

	// R+=c*M*a
	// segment gives the chunk of vector starting at the first argument, and going for as many elements as the second argument...
	// in this case, segment gets the 3vector starting at the 0th DOF's offset (ie 0)
	//R.segment(loadable->GetSubBlockOffset(0), 3) += c * (this->mass * (a_x + chrono::Vcross(a_w, this->c_m))).eigen();
	// in this case, segment gets the 3vector starting at the 0th DOF's + 3 offset (ie 3)
	//R.segment(loadable->GetSubBlockOffset(0) + 3, 3) += c * (this->mass * chrono::Vcross(this->c_m, a_x) + this->I * a_w).eigen();
	// since R is a vector, we can probably just do R += C*M*a with no need to separate w into a_x and a_w above
	R += c * jacobians->M * w;
}