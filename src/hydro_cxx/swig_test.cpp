#include "swig_test.h"
#include <iostream>

void test() {
	std::cout << "Hello World!" << std::endl;
}

// =============================================================================
// HydroInputs Class Definitions
// =============================================================================

/*******************************************************************************
* HydroInputs constructor
*******************************************************************************/
HydroInputs::HydroInputs() {
	// TODO: switch depending on wave option (regular, regularCIC, irregular, noWaveCIC) enum?
	mode = NONE;
	regular_wave_amplitude = 0;
	//excitation_force_phase.resize(6,0);
	//excitation_force_mag.resize(6,0);
}

/*******************************************************************************
* HydroInputs constructor
*******************************************************************************/
HydroInputs::HydroInputs(const HydroInputs& old) {
	*this = old;
}

/*******************************************************************************
* HydroInputs constructor
*******************************************************************************/
HydroInputs& HydroInputs::operator = (const HydroInputs& rhs) {
	freq_index_des = rhs.freq_index_des;
	regular_wave_amplitude = rhs.regular_wave_amplitude;
	regular_wave_omega = rhs.regular_wave_omega;
	wave_omega_delta = rhs.wave_omega_delta;
	excitation_force_mag = rhs.excitation_force_mag;
	excitation_force_phase = rhs.excitation_force_phase;
	mode = rhs.mode;
	return *this;
}
/*******************************************************************************
* HydroInputs test member function for swig wrapping
*******************************************************************************/
void HydroInputs::test2() {
	std::cout << "member function works" << std::endl;
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
TestHydro::TestHydro(const std::vector<std::shared_ptr<ChBody>>& user_bodies, std::string h5_file_name, HydroInputs& user_hydro_inputs) : TestHydro() {
	bodies = user_bodies; // 0 indexed
	num_bodies = bodies.size();
	for (int b = 0; b < num_bodies; b++) {
		file_info.emplace_back(h5_file_name, bodies[b]->GetNameString()); // set up vector of file infos for each body
	}
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

	// added mass info
	my_loadcontainer = chrono_types::make_shared<chrono::ChLoadContainer>();
	my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(file_info, bodies);
	bodies[0]->GetSystem()->Add(my_loadcontainer);
	my_loadcontainer->Add(my_loadbodyinertia);

	// set up hydro inputs stuff 
	hydro_inputs = user_hydro_inputs;
	WaveSetUp();
}

void TestHydro::WaveSetUp() {
	int total_dofs = 6 * num_bodies;
	switch (hydro_inputs.mode) {
	case NONE:

		break;
	case REGULAR:
		hydro_inputs.excitation_force_mag.resize(total_dofs, 0);
		hydro_inputs.excitation_force_phase.resize(total_dofs, 0);
		force_excitation_freq.resize(total_dofs, 0.0);
		hydro_inputs.wave_omega_delta = file_info[0].GetOmegaDelta();
		hydro_inputs.freq_index_des = (hydro_inputs.regular_wave_omega / hydro_inputs.wave_omega_delta) - 1;
		for (int b = 0; b < num_bodies; b++) {
			for (int rowEx = 0; rowEx < 6; rowEx++) {
				int body_offset = 6 * b;
				hydro_inputs.excitation_force_mag[body_offset + rowEx] = file_info[b].GetExcitationMagInterp(rowEx, 0, hydro_inputs.freq_index_des);
				hydro_inputs.excitation_force_phase[body_offset + rowEx] = file_info[b].GetExcitationPhaseInterp(rowEx, 0, hydro_inputs.freq_index_des);
			}
		}
		break;
		// add case: IRREGULAR here TODO
	}
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
	position.resize(total_dofs, 0);
	displacement.resize(total_dofs, 0);

	// orientation initialized to system current pos/rot vectors
	for (int b = 0; b < num_bodies; b++) {
		for (int i = 0; i < 3; i++) {
			unsigned b_offset = 6 * b;
			if (b_offset + i + 3 > total_dofs || b_offset + i < 0) {
				std::cout << "temp index in hydrostatic force is bad " << std::endl;
			}
			//double linearPosition = bodies[b]->GetPos().eigen()[i];
			//double rotationalPosition = bodies[b]->GetRot().Q_to_Euler123()[i];
			position[i + b_offset] = bodies[b]->GetPos().eigen()[i];
			//position[i + 3 + b_offset] = rotationalPosition;
		}
	}

	// make displacement vector for system
	for (int i = 0; i < total_dofs; i++) {
		displacement[i] = position[i] - equilibrium[i];
	}

	// re invent matrix vector multiplication
	for (int b = 0; b < num_bodies; b++) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				//			//KHSOut << i << " " << j << " " << file_info[b].GetHydrostaticStiffness(i, j) << "\n";
				force_hydrostatic[i + (6 * b)] -= ((file_info[b].GetHydrostaticStiffness(i, j)) * displacement[i + (6 * b)]);
			}
		}
	}

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
			setVelHistory(bodies[b - 1]->GetPos_dt()[i],
				(((size + offset_rirf) % size) + size) % size, b, i);
			setVelHistory(bodies[b - 1]->GetWvel_par()[i],
				(((size + offset_rirf) % size) + size) % size, b, i + 3);
		}
	}
	int vi;
	//#pragma omp parallel for
	if (convTrapz == true) {
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
					TMP_S(row, st) = TIMESERIES(row, col, st);
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
* TestHydro::ComputeForceExcitationRegularFreq()
* computes the 6N dimensional excitation force
*******************************************************************************/
std::vector<double> TestHydro::ComputeForceExcitationRegularFreq() {
	for (int b = 0; b < num_bodies; b++) {
		int body_offset = 6 * b;
		for (int rowEx = 0; rowEx < 6; rowEx++) {
			force_excitation_freq[body_offset + rowEx] = hydro_inputs.excitation_force_mag[body_offset + rowEx]
				* hydro_inputs.regular_wave_amplitude * cos(hydro_inputs.regular_wave_omega * bodies[0]->GetChTime()
					+ hydro_inputs.excitation_force_phase[rowEx]);
		}
	}
	return force_excitation_freq;
}

/*******************************************************************************
* TestHydro::TODO
* ****IMPORTANT b is 1 indexed here
* computes all forces or returns saved force if it's already been calculated this
* timestep
* b: body number, it is 1 indexed here since it comes from ForcFunc6d
* i: index or Degree of Freedom (dof) ranges [0,...5]
* calls computeForce type functions
*******************************************************************************/
double TestHydro::coordinateFunc(int b, int i) {
	unsigned body_num_offset = 6 * (b - 1);  // b_num from ForceFunc6d is 1 indexed, TODO: make all b_num 0 indexed
	unsigned total_dofs = 6 * num_bodies;
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

	// reset forces to 0
	std::fill(force_hydrostatic.begin(), force_hydrostatic.end(), 0.0);
	std::fill(force_radiation_damping.begin(), force_radiation_damping.end(), 0);
	std::fill(force_excitation_freq.begin(), force_excitation_freq.end(), 0);

	//call compute forces
	convTrapz = true; // use trapeziodal rule or assume fixed dt.

	ComputeForceHydrostatics();
	ComputeForceRadiationDampingConv();

	// sum all forces element by element
	for (int j = 0; j < total_dofs; j++) {
		total_force[j] = force_hydrostatic[j] - force_radiation_damping[j];
	}

	if (hydro_inputs.mode == REGULAR) {
		ComputeForceExcitationRegularFreq();
		for (int j = 0; j < total_dofs; j++) {
			total_force[j] += force_excitation_freq[j];
		}
	}
	if (body_num_offset + i < 0 || body_num_offset >= total_dofs) {
		std::cout << "total force accessing out of bounds" << std::endl;
	}

	return total_force[body_num_offset + i];
}
