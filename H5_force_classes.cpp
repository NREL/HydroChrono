#include "H5_force_classes.h"
// =============================================================================

/*******************************************************************************
* BodyFileInfo::read_data()
* private member function called from constructor 
* reads h5 file data and stores it in member variables for use with other
* classes and forces
*******************************************************************************/
void BodyFileInfo::read_data() {
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
	double *temp; // TODO change to dynamic memory and delete?
	temp = new double[dims[0] * dims[1]];
	// read file info into data_out, a 2d array
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// turn the 2d array into a ChMatrix (Eigen dynamic matrix)
	lin_matrix.resize(6, 6); // TODO use dims[0] and dims [1]?
	for (int i = 0; i < dims[0]; i++) {
		for (int j = 0; j < dims[1]; j++) {
			lin_matrix(i, j) = temp[i * dims[1] + j];
		}
	}
	lin_matrix *= 10000; // Units are off? scale by 10000? dividing by water density?
	dataset.close();
	delete [] temp;

	// repeat same steps from above to get the cb and cg...reusing some of the previous arrays etc
	data_name = bodyNum + "/properties/cb";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	for (int i = 0; i < 3; i++) {
		cb[i] = temp[i];
	}
	dataset.close();
	//delete [] temp;

	// repeat finally for cg
	data_name = bodyNum + "/properties/cg";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	//temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	for (int i = 0; i < 3; i++) {
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
	// put into equil chvector
	disp_vol = temp[0];
	dataset.close();
	//delete[] temp;

	// read rho
	data_name = "simulation_parameters/rho";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	//temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	rho = temp[0];
	dataset.close();
	//delete[] temp;

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

	// read K
	data_name = bodyNum + "/hydro_coeffs/radiation_damping/impulse_response_fun/K";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	hsize_t dims3[3];    // dataset dimensions
	int rank3 = filespace.getSimpleExtentDims(dims3);
	// read file into data_out 2d array
	H5::DataSpace mspace3(rank3, dims3);
	temp = new double[dims3[0]*dims3[1]*dims3[2]]; // TODO change to dynamic memory and delete?
	// read file info into data_out, a 2d array
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace3, filespace);
	// turn the 2d array into a ChMatrix (Eigen dynamic matrix)
	for (int k = 0; k < dims3[2]; k++) {
		K_matrix[k].resize(6, 6);
		for (int i = 0; i < dims3[0]; i++) {
			for (int j = 0; j < dims3[1]; j++) {
				K_matrix[k](i, j) = temp[ k + dims3[2]*( i * dims3[1] + j ) ]; 
				// h5 file stores first elm of each matrix first, then second elem of each matrix next, row major
			}
		}
	}
	dataset.close();
	delete[] temp;

	data_name = bodyNum + "/hydro_coeffs/radiation_damping/impulse_response_fun/t";
	dataset = sphereFile.openDataSet(data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	temp = new double[dims[0] * dims[1]];
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	for (int i = 0; i < dims[1]; i++) {
		cg[i] = temp[i];
	}
	dataset.close();
	delete[] temp;

	sphereFile.close();
}

BodyFileInfo::BodyFileInfo() {}

/*******************************************************************************
* BodyFileInfo constructor
* requires file name (in absolute file name or referenced from executable location)
* and body name (name of body's section in H5 file, ie "body1" etc)
* each body in system should have its own BodyFileInfo object
* calls read_data()
*******************************************************************************/
BodyFileInfo::BodyFileInfo(std::string file, std::string bodyName) {
	h5_file_name = file;
	bodyNum = bodyName;
	read_data();
}

/*******************************************************************************
* BodyFileInfo::get_lin_matrix()
* returns the linear restoring stiffness matrix 
*******************************************************************************/
ChMatrixDynamic<double> BodyFileInfo::get_lin_matrix() const {
	return lin_matrix;
}

/*******************************************************************************
* BodyFileInfo::get_equil_cg()
* returns cg, center of gravity of object's body
*******************************************************************************/
ChVector<> BodyFileInfo::get_equil_cg() const {
	return cg;
}

/*******************************************************************************
* BodyFileInfo::get_equil_cb()
* returns cb, the rotational equilibrium
*******************************************************************************/
ChVector<> BodyFileInfo::get_equil_cb() const {
	return cb;
}

/*******************************************************************************
* BodyFileInfo::get_rho()
* returns the density of water, rho (kg/m^3 usually)
*******************************************************************************/
double BodyFileInfo::get_rho() const {
	return rho;
}

/*******************************************************************************
* BodyFileInfo::get_g()
* returns g, gravitational acceleration, m/s^2
*******************************************************************************/
double BodyFileInfo::get_g() const {
	return g;
}

/*******************************************************************************
* BodyFileInfo::get_disp_vol() 
* returns displaced volume when body at equilibrium, m^3
*******************************************************************************/
double BodyFileInfo::get_disp_vol() const {
	return disp_vol;
}

/*******************************************************************************
* BodyFileInfo::get_impulse_resp_matrix()
* returns impulse response function matrix K for step i (TODO how to get i)
*******************************************************************************/
ChMatrixDynamic<double> BodyFileInfo::get_impulse_resp_matrix(int i) const {
	return K_matrix[i];
}

// =============================================================================

/*******************************************************************************
* ForceTorqueFunc constructor
* sets pointer to LinRestorForce member object and index for which component
* this ForceTorqueFunc object refers to
*******************************************************************************/
ForceTorqueFunc::ForceTorqueFunc(LinRestorForce* b, int i) : base(b), index(i) {	}

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
void ForceTorqueFunc::SetBase(LinRestorForce* b) {
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

/*******************************************************************************
* LinRestorForce default constructor
* initializes array of ForceTorqueFunc objects and pointers to each force/torque
*******************************************************************************/
LinRestorForce::LinRestorForce() : forces{ {this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5} } {
	for (unsigned i = 0; i < 6; i++) {
		force_ptrs[i] = std::shared_ptr<ForceTorqueFunc>(forces + i, [](ForceTorqueFunc*) {});
		// sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is 
		// automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
		// them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the 
		// default deletion logic to do nothing
		// Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
	}
}

/*******************************************************************************
* LinRestorForce constructor
* calls default constructor and initializes linear restoring stiffness force info
* from BodyFileInfo
* also initializes ChBody that this force will be applied to
*******************************************************************************/
LinRestorForce::LinRestorForce(BodyFileInfo& lin, std::shared_ptr<ChBody> object) : LinRestorForce() {
	//TODO check lin is initialized
	bobber = object;
	fileInfo = lin;
	equil << fileInfo.get_equil_cg().eigen(), fileInfo.get_equil_cb().eigen();
}

/*******************************************************************************
* LinRestorForce::Get_p()
* calculates the matrix multiplication each time step for linear restoring stiffness
* f = [linear restoring stiffness matrix] [displacement vector]
*******************************************************************************/
ChVectorN<double, 6> LinRestorForce::Get_p() const {
	//TODO check this is correct?
	ChVectorN<double, 6> temp;
	temp << bobber->GetPos().eigen(), bobber->GetRot().Q_to_Euler123().eigen();
	temp = equil - temp;
	temp = fileInfo.get_lin_matrix() * temp;
	return temp;
}

/*******************************************************************************
* LinRestorForce::coordinateFunc
* if index is in [0,6] the corresponding vector component of the force vector
* is returned
* otherwise a warning is printed and the force is interpreted to be 0
*******************************************************************************/
double LinRestorForce::coordinateFunc(int i) {
	if (i >= 0 && i < 6) {
		return Get_p()[i];
	}
	else {
		std::cout << "wrong index" << std::endl;
		return 0;
	}
}

/*******************************************************************************
* LinRestorForce::SetForce
* used to initialize components of force (external ChForce pointer)
*******************************************************************************/
void LinRestorForce::SetForce(std::shared_ptr<ChForce> force) {
	force->SetF_x(force_ptrs[0]);
	force->SetF_y(force_ptrs[1]);
	force->SetF_z(force_ptrs[2]);
}

/*******************************************************************************
* LinRestorForce::SetTorque
* used to initialize components of torque (external ChForce pointer with TORQUE flag set)
*******************************************************************************/
void LinRestorForce::SetTorque(std::shared_ptr<ChForce> torque) {
	torque->SetF_x(force_ptrs[3]);
	torque->SetF_y(force_ptrs[4]);
	torque->SetF_z(force_ptrs[5]);
}

// =============================================================================

BuoyancyForce::BuoyancyForce(BodyFileInfo& file) {
	fileInfo = file;
	// get value from file
	bf = fileInfo.get_rho() * fileInfo.get_g() * fileInfo.get_disp_vol();
	// set function to y = bf
	fc.Set_yconst(bf);
	// set pointer to function y=bf
	fc_ptr = std::shared_ptr<ChFunction_Const>(&fc, [](ChFunction_Const*) {} );
	// set force as function pointer in +z direction
	force.SetF_z(fc_ptr);
	// have force_ptr point to force
	force_ptr = std::shared_ptr<ChForce>(&force, [](ChForce*) {});
	// buoyancy force should be [0 0 671980]
}

/*******************************************************************************
*******************************************************************************/
std::shared_ptr<ChForce> BuoyancyForce::getForce_ptr() {
	return force_ptr;
}

/*******************************************************************************
*******************************************************************************/
IRF_func::IRF_func(ImpulseResponseForce* b, int i) : base(b), index(i) { }

/*******************************************************************************
*******************************************************************************/
IRF_func* IRF_func::Clone() const {
	return new IRF_func(*this);
}

/*******************************************************************************
*******************************************************************************/
double 	IRF_func::Get_y(double x) const {
	return base->coordinateFunc(index);
}

/*******************************************************************************
*******************************************************************************/
void IRF_func::SetBase(ImpulseResponseForce* b) {
	base = b;
}

/*******************************************************************************
*******************************************************************************/
void IRF_func::SetIndex(int i) { 
	index = i;
}

/*******************************************************************************
*******************************************************************************/
ImpulseResponseForce::ImpulseResponseForce() : forces{ {this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5} } {
	for (unsigned i = 0; i < 6; i++) {
		force_ptrs[i] = std::shared_ptr<IRF_func>(forces + i, [](IRF_func*) {});
		// sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is 
		// automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
		// them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the 
		// default deletion logic to do nothing
		// Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
	}
}

/*******************************************************************************
*******************************************************************************/
ImpulseResponseForce::ImpulseResponseForce(BodyFileInfo& file, std::shared_ptr<ChBody> object) : ImpulseResponseForce() {
	body = object;
	fileInfo = file;
	// initialize other things from file here
}
/*******************************************************************************
*******************************************************************************/
ChVectorN<double, 6> ImpulseResponseForce::Get_p() const {
	// This function does all the matrix multiplication etc for IRF stuff each timestep
	// should return the combined force/torque 6dof vector to apply to body
	ChVectorN<double, 6> temp;
	for (int i = 0; i < 6; i++) {
		temp[i] = 0;
	}
	return temp;
}

/*******************************************************************************
*******************************************************************************/
double ImpulseResponseForce::coordinateFunc(int i) {
	if (i >= 0 && i < 6) {
		return Get_p()[i];
	}
	else {
		std::cout << "wrong index" << std::endl;
		return 0;
	}
}

/*******************************************************************************
*******************************************************************************/
void ImpulseResponseForce::SetForce(std::shared_ptr<ChForce> force) {
	force->SetF_x(force_ptrs[0]);
	force->SetF_y(force_ptrs[1]);
	force->SetF_z(force_ptrs[2]);
}

/*******************************************************************************
*******************************************************************************/
void ImpulseResponseForce::SetTorque(std::shared_ptr<ChForce> torque) {
	torque->SetF_x(force_ptrs[3]);
	torque->SetF_y(force_ptrs[4]);
	torque->SetF_z(force_ptrs[5]);
}