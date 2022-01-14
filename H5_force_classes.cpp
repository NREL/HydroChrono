#include "H5_force_classes.h"
// =============================================================================

void BodyFileInfo::read_data() {
	// testing HDF5 compatibility
	// open file with read only access
	H5::H5File sphereFile(h5_file_name, H5F_ACC_RDONLY); // "../../test_for_chrono/sphere.h5"

	//
	// Read linear restoring stiffness file info into matrices
	// 
	std::string lin_matrix_data_name = bodyNum + "/hydro_coeffs/linear_restoring_stiffness"; // "body1/hydro_coeffs/linear_restoring_stiffness"
	H5::DataSet dataset = sphereFile.openDataSet(lin_matrix_data_name);
	// Get filespace for rank and dimension
	H5::DataSpace filespace = dataset.getSpace();
	// Get number of dimensions in the file dataspace
	// Get and print the dimension sizes of the file dataspace
	hsize_t dims[2];    // dataset dimensions
	int rank = filespace.getSimpleExtentDims(dims);
	// read file into data_out 2d array
	H5::DataSpace mspace1(rank, dims);
	double temp[36]; // TODO change to dynamic memory and delete?
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

	// repeat same steps from above to get the cb and cg...reusing some of the previous arrays etc
	std::string cb_data_name = bodyNum + "/properties/cb";
	dataset = sphereFile.openDataSet(cb_data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	for (int i = 0; i < 3; i++) {
		cb[i] = temp[i];
	}
	dataset.close();

	// repeat finally for cg
	std::string cg_data_name = bodyNum + "/properties/cg";
	dataset = sphereFile.openDataSet(cg_data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	for (int i = 0; i < 3; i++) {
		cg[i] = temp[i];
	}
	dataset.close();

	// keep reading other things
	
	// read displaced volume for buoyancy force
	std::string disp_vol_data_name = bodyNum + "/properties/disp_vol";
	dataset = sphereFile.openDataSet(disp_vol_data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	// put into equil chvector
	disp_vol = temp[0];
	dataset.close();

	// read rho
	std::string rho_data_name = "simulation_parameters/rho";
	dataset = sphereFile.openDataSet(rho_data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	rho = temp[0];
	dataset.close();

	// read g
	std::string g_data_name = "simulation_parameters/g";
	dataset = sphereFile.openDataSet(g_data_name);
	filespace = dataset.getSpace();
	rank = filespace.getSimpleExtentDims(dims);
	mspace1 = H5::DataSpace(rank, dims);
	dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
	g = temp[0];
	dataset.close();

	sphereFile.close();
}

BodyFileInfo::BodyFileInfo() {}

BodyFileInfo::BodyFileInfo(std::string file, std::string bodyName) {
	h5_file_name = file;
	bodyNum = bodyName;
	read_data();
}

ChMatrixDynamic<double> BodyFileInfo::get_lin_matrix() const {
	return lin_matrix;
}

ChVector<> BodyFileInfo::get_equil_cg() const {
	return cg;
}

ChVector<> BodyFileInfo::get_equil_cb() const {
	return cb;
}

double BodyFileInfo::get_rho() const {
	return rho;
}
double BodyFileInfo::get_g() const {
	return g;
}
double BodyFileInfo::get_disp_vol() const {
	return disp_vol;
}

// =============================================================================

ForceTorqueFunc::ForceTorqueFunc(LinRestorForce* b, int i) : base(b), index(i) {	}

ForceTorqueFunc* ForceTorqueFunc::Clone() const  { 
	return new ForceTorqueFunc(*this); 
}

double ForceTorqueFunc::Get_y(double x) const {
	return base->coordinateFunc(index);
}

void ForceTorqueFunc::SetBase(LinRestorForce* b) {
	base = b;
}

void ForceTorqueFunc::SetIndex(int i) {
	index = i;
}

// =============================================================================

LinRestorForce::LinRestorForce() : forces{ {this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5} } {
	for (unsigned i = 0; i < 6; i++) {
		force_ptrs[i] = std::shared_ptr<ForceTorqueFunc>(forces + i, [](ForceTorqueFunc*) {});
	}
}

LinRestorForce::LinRestorForce(BodyFileInfo& lin, std::shared_ptr<ChBody> object) : LinRestorForce() {
	//TODO check lin is initialized
	bobber = object;
	fileInfo = lin;
	equil << fileInfo.get_equil_cg().eigen(), fileInfo.get_equil_cb().eigen();
}

ChVectorN<double, 6> LinRestorForce::Get_p() const {
	//TODO check this is correct?
	ChVectorN<double, 6> temp;
	temp << bobber->GetPos().eigen(), bobber->GetRot().Q_to_Euler123().eigen();
	temp = equil - temp;
	temp = fileInfo.get_lin_matrix() * temp;
	return temp;
}

double LinRestorForce::coordinateFunc(int i) {
	if (i >= 0 && i < 6) {
		return Get_p()[i];
	}
	else {
		std::cout << "wrong index" << std::endl;
		return 0;
	}
}

void LinRestorForce::SetForce(std::shared_ptr<ChForce> force) {
	force->SetF_x(force_ptrs[0]);
	force->SetF_y(force_ptrs[1]);
	force->SetF_z(force_ptrs[2]);
}

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

std::shared_ptr<ChForce> BuoyancyForce::getForce_ptr() {
	return force_ptr;
}

