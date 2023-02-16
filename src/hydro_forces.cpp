#include "hydroc/hydro_forces.h"

#include <algorithm>
#include <numeric>  // std::accumulate

#include <hydroc/h5fileinfo.h>
#include <hydroc/chloadaddedmass.h>

// =============================================================================
// HydroInputs Class Definitions
// =============================================================================

/*******************************************************************************
 * HydroInputs constructor
 *******************************************************************************/
HydroInputs::HydroInputs() {
    // TODO: switch depending on wave option (regular, regularCIC, irregular, noWaveCIC) enum?
    mode = WaveMode::noWaveCIC;
    regular_wave_amplitude = 0.0;
    excitation_force_phase.resize(6, 0);
    excitation_force_mag.resize(6, 0);
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
ComponentFunc::ComponentFunc(ForceFunc6d* b, int i) : base(b), index(i) {}

/*******************************************************************************
 * ComponentFunc::Clone()
 * required override function since ComponentFunc inherits from ChFunction
 *******************************************************************************/
ComponentFunc* ComponentFunc::Clone() const {
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
ForceFunc6d::ForceFunc6d() : forces{{this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5}} {
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
    chrono_force->SetAlign(ChForce::AlignmentFrame::WORLD_DIR);
    chrono_torque->SetAlign(ChForce::AlignmentFrame::WORLD_DIR);
    chrono_force->SetNameString("hydroforce");
    chrono_torque->SetNameString("hydrotorque");
}

/*******************************************************************************
 * ForceFunc6d constructor
 * calls default constructor and initializes hydro force info
 * from H5FileInfo
 * also initializes ChBody that this force will be applied to
 *******************************************************************************/
ForceFunc6d::ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* user_all_forces) : ForceFunc6d() {
    body = object;
    std::string temp = body->GetNameString();  // remove "body" from "bodyN", convert N to int, get body num
    b_num = stoi(temp.erase(0, 4));            // 1 indexed
    all_hydro_forces = user_all_forces;        // TODO switch to smart pointers? does this use = ?
    if (all_hydro_forces == NULL) {
        std::cout << "all hydro forces null " << std::endl;
    }
    SetForce();
    SetTorque();
    ApplyForceAndTorqueToBody();
}

/*******************************************************************************
 * ForceFunc6d copy constructor
 * copy constructor should check to see if this force has been added to this body yet
 * if not, it should add it, if so it shouldnt add the force a second time
 * TODO
 *******************************************************************************/
ForceFunc6d::ForceFunc6d(const ForceFunc6d& old)
    : forces{{this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5}} {
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
}

/*******************************************************************************
 * ForceFunc6d::ApplyForceAndTorqueToBody()
 * adds this force to the body's list of applied forces
 * Warning: everytime this is called, a force is applied to the body so be careful not to
 * duplicate forces on accident
 * TODO: make this function less risky, there shouldn't be a chance to apply the same force
 * multiple times
 *******************************************************************************/
void ForceFunc6d::ApplyForceAndTorqueToBody() {
    body->AddForce(chrono_force);
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
    prev_time = -1;
    offset_rirf = 0;
    num_bodies = 0;
}

#include <chrono/physics/ChLoad.h>
#include <vector>
#include <memory>

/*******************************************************************************
 * TestHydro::TestHydro(user_bodies, h5_file_name, user_hydro_inputs)
 * main constructor for TestHydro class, sets up vector of bodies, h5 file info,
 * and hydro inputs
 * also initializes many persistent variables for force calculations
 * calls default constructor
 *******************************************************************************/
TestHydro::TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
                     std::string h5_file_name,
                     HydroInputs& user_hydro_inputs)
    : TestHydro() {
    bodies = user_bodies;  // 0 indexed
    num_bodies = bodies.size();

    for (int b = 0; b < num_bodies; b++) {
        file_info.emplace_back(h5_file_name, bodies[b]->GetNameString());  // set up vector of file infos for each body
    }

    // set up time vector (should be the same for each body, so just use the first always)
    rirf_time_vector = file_info[0].GetRIRFTimeVector();
    rirf_timestep = rirf_time_vector[1] - rirf_time_vector[0];  // TODO is this the same for all bodies?
    // simplify 6* num_bodies to be the system's total number of dofs, makes expressions later easier to read
    unsigned total_dofs = 6 * num_bodies;
    // resize and initialize velocity history vector to all zeros
    velocity_history.resize(file_info[0].GetRIRFDims(2) * total_dofs, 0.0);  // resize and fill with 0s
    // resize and initialize all persistent forces to all 0s
    force_hydrostatic.resize(total_dofs, 0.0);
    force_radiation_damping.resize(total_dofs, 0.0);
    total_force.resize(total_dofs, 0.0);
    // set up equilibrium for entire system (each body has position and rotation equilibria 3 indicies apart)
    equilibrium.resize(total_dofs, 0.0);
    cb_minus_cg.resize(3 * num_bodies, 0.0);  // cb-cg has 3 components for each body
    for (int b = 0; b < num_bodies; b++) {
        for (int i = 0; i < 3; i++) {
            unsigned equilibrium_idx = i + 6 * b;
            unsigned c_idx = i + 3 * b;
            equilibrium[equilibrium_idx] = file_info[b].cg[i];  // positional equilib is cg, leave rotational bit 0
            cb_minus_cg[c_idx] = file_info[b].cb[i] - file_info[b].cg[i];
        }
    }

    for (int b = 0; b < num_bodies; b++) {
        force_per_body.emplace_back(bodies[b], this);
    }

    // added mass info
    my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();

    /// TODO Check if local vector is really copied into constructor of ChLoadAddedMass
    /// else it could be a memory fault
    std::vector<std::shared_ptr<ChLoadable>> loadables(bodies.size());
    for (auto i = 0; i < bodies.size(); ++i) {
        loadables[i] = bodies[i];
    }

    my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(file_info, loadables);
    bodies[0]->GetSystem()->Add(my_loadcontainer);
    my_loadcontainer->Add(my_loadbodyinertia);

    // set up hydro inputs stuff
    hydro_inputs = user_hydro_inputs;
    WaveSetUp();
}

void TestHydro::WaveSetUp() {
    int total_dofs = 6 * num_bodies;
    switch (hydro_inputs.mode) {
        case WaveMode::noWaveCIC:
            break;
        case WaveMode::regular:
            hydro_inputs.excitation_force_mag.resize(total_dofs, 0.0);
            hydro_inputs.excitation_force_phase.resize(total_dofs, 0.0);
            force_excitation_freq.resize(total_dofs, 0.0);
            hydro_inputs.wave_omega_delta = file_info[0].GetOmegaDelta();
            hydro_inputs.freq_index_des = (hydro_inputs.regular_wave_omega / hydro_inputs.wave_omega_delta) - 1;
            for (int b = 0; b < num_bodies; b++) {
                for (int rowEx = 0; rowEx < 6; rowEx++) {
                    int body_offset = 6 * b;
                    hydro_inputs.excitation_force_mag[body_offset + rowEx] =
                        file_info[b].GetExcitationMagInterp(rowEx, 0, hydro_inputs.freq_index_des);
                    hydro_inputs.excitation_force_phase[body_offset + rowEx] =
                        file_info[b].GetExcitationPhaseInterp(rowEx, 0, hydro_inputs.freq_index_des);
                }
            }
            break;
            // add case: IRREGULAR here TODO
    }
}

/*******************************************************************************
 * TestHydro::getVelHistoryVal(int step, int c) const
 * finds and returns the component of velocity history for given step and dof (c - column)
 * step: [0,1,...,1000] (timesteps from h5 file, one velocity per step
 * c: [0,..,num_bodies-1,...,numbodies*6-1] (in order of bodies, iterates over
 *    dof for each body...3 bodies c would be [0,1,...,17])
 *******************************************************************************/
double TestHydro::getVelHistoryVal(int step, int c) const {
    if (step < 0 || step >= file_info[0].GetRIRFDims(2) || c < 0 || c >= num_bodies * 6) {
        std::cout << "wrong vel history index " << std::endl;
        return 0;
    }
    int index = c % 6;
    int b = c / 6;  // 0 indexed
    if (index + 6 * b + 6 * num_bodies * step >= num_bodies * 6 * file_info[0].GetRIRFDims(2) ||
        index + 6 * b + 6 * num_bodies * step < 0) {
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
    if (index + 6 * (b_num - 1) + 6 * num_bodies * step < 0 ||
        index + 6 * (b_num - 1) + 6 * num_bodies * step >= num_bodies * 6 * file_info[0].GetRIRFDims(2)) {
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
    position.resize(total_dofs, 0.0);
    displacement.resize(total_dofs, 0.0);

    // orientation initialized to system current pos/rot vectors
    for (int b = 0; b < num_bodies; b++) {
        for (int i = 0; i < 3; i++) {
            unsigned b_offset = 6 * b;
            if (b_offset + i + 3 > total_dofs || b_offset + i < 0) {
                std::cout << "temp index in hydrostatic force is bad " << std::endl;
            }
            position[i + b_offset] = bodies[b]->GetPos()[i];
            position[i + 3 + b_offset] =
                bodies[b]->GetRot().Q_to_Euler123()[i];  // TODO check if we want this, or GetWvel_par()?
        }
    }

    // make displacement vector for system
    for (int i = 0; i < total_dofs; i++) {
        displacement[i] = equilibrium[i] - position[i];
    }

    // re invent matrix vector multiplication
    for (int b = 0; b < num_bodies; b++) {
        unsigned b_offset = 6 * b;
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                double t = (file_info[b].GetHydrostaticStiffness(i, j));
                // if (i == 4 && j == 4 && t < 0) {
                //	t = -t;
                // }
                force_hydrostatic[i + b_offset] += t * displacement[j + b_offset];
            }
        }
    }

    // now handle buoyancy force....
    assert(num_bodies > 0);
    double* buoyancy = new double[num_bodies];  // this sets up an array for buoyancy for each body
    // add heave buoyancy for each body, and add rxb=(cb-cg)x(0,0,buoyancy) for the moment due to buoyancy for each body
    // (simplified)
    for (int b = 0; b < num_bodies; b++) {
        buoyancy[b] = file_info[b].rho * -(bodies[b]->GetSystem()->Get_G_acc()).z() *
                      file_info[b].disp_vol;  // buoyancy = rho*g*Vdisp
        unsigned b_offset =
            6 * b;  // force_hydrostatic has 6 elements for each body so to skip to the next body we move 6 spaces
        unsigned r_offset =
            3 * b;  // cb_minus_cg has 3 elements for each body so to skip to the next body we move 3 spaces
        force_hydrostatic[b_offset + 2] += buoyancy[b];  // add heave buoyancy
        // now for moments due to buoyancy (simplified)
        // for torque about x (index 3) per body, add b * r_y
        force_hydrostatic[b_offset + 3] += buoyancy[b] * cb_minus_cg[1 + r_offset];
        // for torque about y (index 4) per body, subtract b*r_x
        force_hydrostatic[b_offset + 4] += -1 * buoyancy[b] * cb_minus_cg[0 + r_offset];
    }
    delete[] buoyancy;
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
    offset_rirf--;  // starts as 0 before timestep change
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
#define TIMESERIES(row, col, step) timeseries[(row)*numCols * size + (col)*size + (step)]
    // TMP_S ends up being a sum over the columns of TIMESERIES (total_dofs aka LDOF
#define TMP_S(row, step) tmp_s[(row)*size + (step)]
    // set last entry as velocity
    for (int i = 0; i < 3; i++) {
        for (int b = 1; b < num_bodies + 1; b++) {  // body index being 1 indexed here is right
            int vi = (((size + offset_rirf) % size) + size) % size;
            setVelHistory(bodies[b - 1]->GetPos_dt()[i], vi, b, i);
            setVelHistory(bodies[b - 1]->GetWvel_par()[i], vi, b, i + 3);
        }
    }
    int vi;
    //#pragma omp parallel for
    if (convTrapz == true) {
        // convolution integral using trapezoidal rule
        for (int row = 0; row < numRows; row++) {  // row goes to 6N
            for (int st = 0; st < size; st++) {
                vi = (((st + offset_rirf) % size) + size) % size;  // vi takes care of circshift function from matLab
                TMP_S(row, st) = 0;
                for (int col = 0; col < numCols; col++) {  // numCols goes to 6N
                    // multiply rirf by velocity history for each step and row (0,...,6N), store product in TIMESERIES
                    TIMESERIES(row, col, st) = GetRIRFval(row, col, st) * getVelHistoryVal(vi, col);
                    // for (int i = 0; i < numCols; i++) {
                    // velOut << getVelHistoryVal(vi, col) << std::endl;
                    // }
                    //  TMP_S is the sum over col (sum the effects of all radiating dofs (LDOF) for each time and motion
                    //  dof)
                    TMP_S(row, st) += TIMESERIES(row, col, st);
                }
                if (st > 0) {
                    // integrate TMP_S
                    force_radiation_damping[row] +=
                        (TMP_S(row, st - 1) + TMP_S(row, st)) / 2.0 * (rirf_time_vector[st] - rirf_time_vector[st - 1]);
                }
            }
        }
    }
    // velOut.close();

    // else { // TODO fix this for force_radiation_damping to go over col not row!
    //	// convolution integral assuming fixed dt
    //	for (int row = 0; row < numRows; row++) {
    //		//#pragma omp parallel for
    //		sumVelHistoryAndRIRF = 0.0;
    //		for (int col = 0; col < numCols; col++) {
    //			for (int st = 0; st < size; st++) {
    //				vi = (((st + offset_rirf) % size) + size) % size; // vi takes care of circshift function from matLab
    //				TIMESERIES(row, col, st) = GetRIRFval(row, col, st) * getVelHistoryVal(vi, col); // col now runs thru all
    //bodies (0->11 for 2 bodies...) 				TMP_S(row, st) = TIMESERIES(row, col, st); 				sumVelHistoryAndRIRF += TMP_S(row, st);
    //			}
    //		}
    //		force_radiation_damping[row] -= sumVelHistoryAndRIRF * rirf_timestep;
    //	}
    // }
    //  Deallocate memory
#undef TIMESERIES
#undef TMP_S
    delete[] timeseries;
    delete[] tmp_s;

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
    if (row < 0 || row >= 6 * num_bodies || col < 0 || col >= 6 * num_bodies || st < 0 ||
        st >= file_info[0].GetRIRFDims(2)) {
        std::cout << "rirfval index bad from testhydro" << std::endl;
        return 0;
    }
    int b = row / 6;  // 0 indexed, which body to get matrix info from
    int c = col % 6;  // which dof across column, 0,..,11 for 2 bodies, 0,...,6N-1 for N
    int r = row % 6;  // which dof 0,..,5 in individual body RIRF matrix
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
            force_excitation_freq[body_offset + rowEx] = hydro_inputs.excitation_force_mag[body_offset + rowEx] *
                                                         hydro_inputs.regular_wave_amplitude *
                                                         cos(hydro_inputs.regular_wave_omega * bodies[0]->GetChTime() +
                                                             hydro_inputs.excitation_force_phase[rowEx]);
        }
    }
    return force_excitation_freq;
}

///*******************************************************************************
//* TestHydro::ComputeForceExcitationRegularFreq()
//* computes the 6N dimensional excitation force
//*******************************************************************************/
// std::vector<double> TestHydro::ComputeForceRegularWaves() {
//	//for (int b = 0; b < num_bodies; b++) {
//	//	int body_offset = 6 * b;
//	//	for (int rowEx = 0; rowEx < 6; rowEx++) {
//	//		if (rowEx == 2) {
//	//			force_excitation_freq[body_offset + rowEx] = hydro_inputs.excitation_force_mag[rowEx]
//	//				* hydro_inputs.regular_wave_amplitude * cos(hydro_inputs.regular_wave_omega * bodies[0]->GetChTime()
//	//					+ excitation_force_phase[rowEx]);
//	//		}
//	//		else {
//	//			force_excitation_freq[rowEx] = 0.0;
//	//		}
//	//	}
//	//}
//	return force_reg_waves;
//}

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
    std::fill(total_force.begin(), total_force.end(), 0.0);
    std::fill(force_hydrostatic.begin(), force_hydrostatic.end(), 0.0);
    std::fill(force_radiation_damping.begin(), force_radiation_damping.end(), 0.0);
    std::fill(force_excitation_freq.begin(), force_excitation_freq.end(), 0.0);

    // call compute forces
    convTrapz = true;  // use trapeziodal rule or assume fixed dt.

    if (hydro_inputs.mode == WaveMode::noWaveCIC) {
        // update required forces:
        ComputeForceHydrostatics();

        std::ofstream hydrostaticsCheck;
        hydrostaticsCheck.open("./results/f3of/debugging/fhs.txt");
        for (int dof = 0; dof < total_dofs; dof++) {
            hydrostaticsCheck << force_hydrostatic[dof] << std::endl;
        }
        hydrostaticsCheck.close();

        ComputeForceRadiationDampingConv();

        std::ofstream radiationCheck;
        radiationCheck.open("./results/f3of/debugging/frad.txt");
        for (int dof = 0; dof < total_dofs; dof++) {
            radiationCheck << force_radiation_damping[dof] << std::endl;
        }
        radiationCheck.close();

        // sum all forces element by element
        for (int j = 0; j < total_dofs; j++) {
            total_force[j] = force_hydrostatic[j] - force_radiation_damping[j];  // force_hydrostatic[j] -
        }

        std::ofstream totalForceCheck;
        totalForceCheck.open("./results/f3of/debugging/ftotal.txt");
        for (int dof = 0; dof < total_dofs; dof++) {
            totalForceCheck << total_force[dof] << std::endl;
        }
        totalForceCheck.close();
    } else if (hydro_inputs.mode == WaveMode::regular) {
        // update required forces:
        ComputeForceHydrostatics();
        ComputeForceRadiationDampingConv();
        ComputeForceExcitationRegularFreq();
        // sum all forces element by element
        for (int j = 0; j < total_dofs; j++) {
            total_force[j] = force_hydrostatic[j] - force_radiation_damping[j] + force_excitation_freq[j];
        }
    }
    if (body_num_offset + i < 0 || body_num_offset >= total_dofs) {
        std::cout << "total force accessing out of bounds" << std::endl;
    }

    return total_force[body_num_offset + i];
}
