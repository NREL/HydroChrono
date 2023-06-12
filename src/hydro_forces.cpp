/*********************************************************************
 * @file  hydro_forces.cpp
 *
 * @brief Implementation of TestHydro main class and helper classes
 * ComponentFunc and ForceFunc6d.
 *********************************************************************/

// TODO minimize include statements, move all to header file hydro_forces.h?
#include "hydroc/hydro_forces.h"
#include <hydroc/chloadaddedmass.h>
#include <hydroc/h5fileinfo.h>
#include <hydroc/wave_types.h>

#include <chrono/physics/ChLoad.h>
#include <unsupported/Eigen/Splines>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>  // std::accumulate
#include <random>
#include <vector>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <vector>

// TODO: move Misc writecontainer type functions to different file
// TODO move WriteContainerToFile generic declaration to a .h file instead of .cpp
/**
 * @brief Prints contents of 1D Container data to given file.
 *
 * @param container the 1D array/vector to write to file
 * @param file_name file to write container to
 */
template <typename Container>
void WriteContainerToFile(const Container& container, const std::string& file_name);

/**
 * @brief Prints contents of std::vector<double> data to given file.
 *
 * @param container std::vector<double> to write to file
 * @param file_name file to write container to
 */
template <>
void WriteContainerToFile<std::vector<double>>(const std::vector<double>& container, const std::string& file_name) {
    std::ofstream output_file(file_name);

    if (!output_file) {
        std::cerr << "Error: Unable to open the file: " << file_name << std::endl;
        return;
    }

    for (const double value : container) {
        output_file << value << std::endl;
    }

    output_file.close();
}

/**
 * @brief Prints contents of Eigen::VectorXd data to given file.
 *
 * @param container Eigen::VectorXd to write to file
 * @param file_name file to write container to
 */
template <>
void WriteContainerToFile<Eigen::VectorXd>(const Eigen::VectorXd& container, const std::string& file_name) {
    std::ofstream output_file(file_name);

    if (!output_file) {
        std::cerr << "Error: Unable to open the file: " << file_name << std::endl;
        return;
    }

    for (int i = 0; i < container.size(); ++i) {
        output_file << container[i] << std::endl;
    }

    output_file.close();
}

// TODO: on move from std::vector to Eigen::VectroXd, remove this function, look up
// Eigen::LinSpaced function as replacement (arguments in slightly different order)
/**
 * @brief Generates a std::vector<double> with evenly spaced points
 *
 * @param start value for the 0th item in the std::vector<double>
 * @param end value for the last item in the std::vector<double>
 * @param num_points defines the size for the std::vector<double>
 *
 * @return std::vector<double> with evenly spaced points from start to end
 */
std::vector<double> Linspace(double start, double end, int num_points) {
    std::vector<double> result(num_points);
    double step = (end - start) / (num_points - 1);

    for (int i = 0; i < num_points; ++i) {
        result[i] = start + i * step;
    }

    return result;
}

// TODO reorder ComponentFunc implementation functions to match the header order of functions
ComponentFunc::ComponentFunc() {
    base  = NULL;
    index = 6;
}

ComponentFunc::ComponentFunc(ForceFunc6d* b, int i) : base(b), index(i) {}

ComponentFunc* ComponentFunc::Clone() const {
    return new ComponentFunc(*this);
}

ComponentFunc::ComponentFunc(const ComponentFunc& old) {
    base  = old.base;
    index = old.index;
}

double ComponentFunc::Get_y(double x) const {
    if (base == NULL) {
        std::cout << "base == Null!" << std::endl;
        return 0;
    }
    return base->coordinateFunc(index);
}

ForceFunc6d::ForceFunc6d() : forces{{this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5}} {
    for (unsigned i = 0; i < 6; i++) {
        force_ptrs[i] = std::shared_ptr<ComponentFunc>(forces + i, [](ComponentFunc*) {});
        // sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is
        // automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
        // them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the
        // default deletion logic to do nothing
        // Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
    }
    chrono_force  = chrono_types::make_shared<ChForce>();
    chrono_torque = chrono_types::make_shared<ChForce>();
    chrono_force->SetAlign(ChForce::AlignmentFrame::WORLD_DIR);
    chrono_torque->SetAlign(ChForce::AlignmentFrame::WORLD_DIR);
    chrono_force->SetNameString("hydroforce");
    chrono_torque->SetNameString("hydrotorque");
}

ForceFunc6d::ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* user_all_forces) : ForceFunc6d() {
    body             = object;
    std::string temp = body->GetNameString();   // remove "body" from "bodyN", convert N to int, get body num
    b_num            = stoi(temp.erase(0, 4));  // 1 indexed TODO: fix b_num starting here to be 0 indexed
    all_hydro_forces = user_all_forces;         // TODO switch to smart pointers? does this use = ?
    if (all_hydro_forces == NULL) {
        std::cout << "all hydro forces null " << std::endl;
    }
    SetForce();
    SetTorque();
    ApplyForceAndTorqueToBody();
}

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
    chrono_force     = old.chrono_force;
    chrono_torque    = old.chrono_torque;
    body             = old.body;
    b_num            = old.b_num;
    all_hydro_forces = old.all_hydro_forces;
    SetForce();
    SetTorque();
}

double ForceFunc6d::coordinateFunc(int i) {
    // b_num is 1 indexed?
    if (i >= 6 || i < 0) {
        std::cout << "wrong index force func 6d" << std::endl;
        return 0;
    }
    return all_hydro_forces->coordinateFunc(
        b_num, i);  // b_num is 1 indexed here!!!!! TODO: change all b_num to be 0 indexed everywhere
}

void ForceFunc6d::SetForce() {
    if (chrono_force == NULL || body == NULL) {
        std::cout << "set force null issue" << std::endl;
    }
    chrono_force->SetF_x(force_ptrs[0]);
    chrono_force->SetF_y(force_ptrs[1]);
    chrono_force->SetF_z(force_ptrs[2]);
}

void ForceFunc6d::SetTorque() {
    if (chrono_torque == NULL || body == NULL) {
        std::cout << "set torque null issue" << std::endl;
    }
    chrono_torque->SetF_x(force_ptrs[3]);
    chrono_torque->SetF_y(force_ptrs[4]);
    chrono_torque->SetF_z(force_ptrs[5]);
    chrono_torque->SetMode(ChForce::ForceType::TORQUE);
}

void ForceFunc6d::ApplyForceAndTorqueToBody() {
    body->AddForce(chrono_force);
    body->AddForce(chrono_torque);
}

// TODO reorder TestHydro function ordering to match header file order
TestHydro::TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
                     std::string h5_file_name,
                     std::shared_ptr<WaveBase> waves)
    : bodies(user_bodies), num_bodies(bodies.size()), file_info(H5FileInfo(h5_file_name, num_bodies).readH5Data()) {
    prev_time   = -1;
    offset_rirf = 0;

    // set up time vector (should be the same for each body, so just use the first always)
    rirf_time_vector = file_info.GetRIRFTimeVector();
    rirf_timestep    = rirf_time_vector[1] - rirf_time_vector[0];  // TODO is this the same for all bodies?

    // simplify 6* num_bodies to be the system's total number of dofs, makes expressions later easier to read
    int total_dofs = 6 * num_bodies;
    // resize and initialize velocity history vector to all zeros
    velocity_history.resize(file_info.GetRIRFDims(2) * total_dofs, 0.0);  // resize and fill with 0s
    // resize and initialize all persistent forces to all 0s
    // TODO rephrase for Eigen::VectorXd eventually
    force_hydrostatic.resize(total_dofs, 0.0);
    force_radiation_damping.resize(total_dofs, 0.0);
    total_force.resize(total_dofs, 0.0);
    // set up equilibrium for entire system (each body has position and rotation equilibria 3 indicies apart)
    equilibrium.resize(total_dofs, 0.0);
    cb_minus_cg.resize(3 * num_bodies, 0.0);  // cb-cg has 3 components for each body
    for (int b = 0; b < num_bodies; b++) {
        for (int i = 0; i < 3; i++) {
            unsigned equilibrium_idx = i + 6 * b;
            unsigned c_idx           = i + 3 * b;
            // positional equilib is cg, leave rotational bit 0
            equilibrium[equilibrium_idx] = file_info.GetCGVector(b)[i];
            cb_minus_cg[c_idx]           = file_info.GetCBVector(b)[i] - file_info.GetCGVector(b)[i];
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

    my_loadbodyinertia =
        chrono_types::make_shared<ChLoadAddedMass>(file_info.GetBodyInfos(), loadables, bodies[0]->GetSystem());
    bodies[0]->GetSystem()->Add(my_loadcontainer);
    my_loadcontainer->Add(my_loadbodyinertia);

    // set up hydro inputs stuff
    // hydro_inputs = user_hydro_inputs;
    // WaveSetUp();
    user_waves = waves;
    AddWaves(user_waves);
}

void TestHydro::AddWaves(std::shared_ptr<WaveBase> waves) {
    user_waves = waves;
    if (user_waves->GetWaveMode() == WaveMode::regular) {
        std::shared_ptr<RegularWave> reg = std::static_pointer_cast<RegularWave>(user_waves);
        reg->AddH5Data(file_info.GetRegularWaveInfos());
    } else if (user_waves->GetWaveMode() == WaveMode::irregular) {
        std::shared_ptr<IrregularWave> irreg = std::static_pointer_cast<IrregularWave>(user_waves);
        irreg->AddH5Data(file_info.GetIrregularWaveInfos(), file_info.GetSimulationInfo());
    }
    user_waves->Initialize();
}

double TestHydro::getVelHistoryVal(int step, int c) const {
    if (step < 0 || step >= file_info.GetRIRFDims(2) || c < 0 || c >= num_bodies * 6) {
        std::cout << "wrong vel history index " << std::endl;
        return 0;
    }
    int index = c % 6;
    int b     = c / 6;  // 0 indexed
    if (index + 6 * b + 6 * num_bodies * step >= num_bodies * 6 * file_info.GetRIRFDims(2) ||
        index + 6 * b + 6 * num_bodies * step < 0) {
        std::cout << "bad vel history math" << std::endl;
        return 0;
    }
    return velocity_history[index + (6 * b) + (6 * num_bodies * step)];
}

double TestHydro::setVelHistory(double val, int step, int b_num, int index) {
    if (step < 0 || step >= file_info.GetRIRFDims(2) || b_num < 1 || b_num > num_bodies || index < 0 || index >= 6) {
        std::cout << "bad set vel history indexing" << std::endl;
        return 0;
    }
    if (index + 6 * (b_num - 1) + 6 * num_bodies * step < 0 ||
        index + 6 * (b_num - 1) + 6 * num_bodies * step >= num_bodies * 6 * file_info.GetRIRFDims(2)) {
        std::cout << "bad set vel history math" << std::endl;
        return 0;
    }
    velocity_history[index + (6 * (b_num - 1)) + (6 * num_bodies * step)] = val;
    return val;
}

std::vector<double> TestHydro::ComputeForceHydrostatics() {
    assert(num_bodies > 0);

    for (int b = 0; b < num_bodies; b++) {
        // initialize variables
        std::shared_ptr<chrono::ChBody> body = bodies[b];
        // H5FileInfo& body_h5file              = file_info[b];
        double rho   = file_info.GetRhoVal();
        int b_offset = 6 * b;
        // force_hydrostatic has 6 elements for each body so to skip to the next body we move 6 spaces
        double* body_force_hydrostatic = &force_hydrostatic[b_offset];
        double* body_equilibrium       = &equilibrium[b_offset];
        double gg                      = body->GetSystem()->Get_G_acc().Length();

        // hydrostatic stiffness due to offset from equilibrium
        chrono::ChVector<> body_position = body->GetPos();
        chrono::ChVector<> body_rotation = body->GetRot().Q_to_Euler123();
        // calculate displacement
        chrono::ChVectorN<double, 6> body_displacement;
        for (int ii = 0; ii < 3; ii++) {
            body_displacement[ii]     = body_position[ii] - body_equilibrium[ii];
            body_displacement[ii + 3] = body_rotation[ii] - body_equilibrium[ii + 3];
        }
        // calculate force
        chrono::ChVectorN<double, 6> force_offset = -gg * rho * file_info.GetLinMatrix(b) * body_displacement;
        // add to force_hydrostatic
        for (int dof = 0; dof < 6; dof++) {
            body_force_hydrostatic[dof] += force_offset[dof];
        }

        // buoyancy at equilibrium
        // TODO: move to prestep (shouldn't be calculated at each time step)
        // translational
        chrono::ChVector<> buoyancy =
            rho * (-body->GetSystem()->Get_G_acc()) * file_info.GetDispVolVal(b);  // buoyancy = rho*g*Vdisp
        body_force_hydrostatic[0] += buoyancy[0];
        body_force_hydrostatic[1] += buoyancy[1];
        body_force_hydrostatic[2] += buoyancy[2];
        // rotational
        int r_offset = 3 * b;
        // cb_minus_cg has 3 elements for each body so to skip to the next body we move 3 spaces
        auto cg2cb =
            chrono::ChVector<double>(cb_minus_cg[0 + r_offset], cb_minus_cg[1 + r_offset], cb_minus_cg[2 + r_offset]);
        chrono::ChVector<> buoyancy2 = cg2cb % buoyancy;
        body_force_hydrostatic[3] += buoyancy2[0];
        body_force_hydrostatic[4] += buoyancy2[1];
        body_force_hydrostatic[5] += buoyancy2[2];
    }
    return force_hydrostatic;
}

std::vector<double> TestHydro::ComputeForceRadiationDampingConv() {
    int size = file_info.GetRIRFDims(2);
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
    double* tmp_s      = new double[numRows * size];
    // define shortcuts for accessing 1D arrays as 3D (or 2D) arrays
    // TIMESERIES is for each row in RIRF, element wise multipy velocity history by RIRF slab
#define TIMESERIES(row, col, step) timeseries[(row * numCols * size) + (col * size) + (step)]
    // TMP_S ends up being a sum over the columns of TIMESERIES (total_dofs aka LDOF
#define TMP_S(row, step) tmp_s[((row)*size) + (step)]
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
    if (convTrapz == true) {  // TODO add public function to TestHydro for users to change the value of convTrapz from
                              // main simulation
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
                    //}
                    // TMP_S is the sum over col (sum the effects of all radiating dofs (LDOF) for each time and motion
                    // dof)
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

    // else { // TODO fix this for force_radiation_damping to go over col not row like above!
    //	// convolution integral assuming fixed dt
    //	for (int row = 0; row < numRows; row++) {
    //		//#pragma omp parallel for
    //		sumVelHistoryAndRIRF = 0.0;
    //		for (int col = 0; col < numCols; col++) {
    //			for (int st = 0; st < size; st++) {
    //				vi = (((st + offset_rirf) % size) + size) % size; // vi takes care of circshift function from matLab
    //				TIMESERIES(row, col, st) = GetRIRFval(row, col, st) * getVelHistoryVal(vi, col); // col now runs
    // thru all
    // bodies (0->11 for 2 bodies...) 				TMP_S(row, st) = TIMESERIES(row, col, st);
    // sumVelHistoryAndRIRF
    // += TMP_S(row, st);
    //			}
    //		}
    //		force_radiation_damping[row] -= sumVelHistoryAndRIRF * rirf_timestep;
    //	}
    //}
    // Deallocate memory
#undef TIMESERIES
#undef TMP_S
    delete[] timeseries;
    delete[] tmp_s;

    return force_radiation_damping;
}

double TestHydro::GetRIRFval(int row, int col, int st) {
    if (row < 0 || row >= 6 * num_bodies || col < 0 || col >= 6 * num_bodies || st < 0 ||
        st >= file_info.GetRIRFDims(2)) {
        std::cout << "rirfval index bad from testhydro" << std::endl;
        return 0;
    }
    int b = row / 6;  // 0 indexed, which body to get matrix info from
    int c = col % 6;  // which dof across column, 0,..,11 for 2 bodies, 0,...,6N-1 for N
    int r = row % 6;  // which dof 0,..,5 in individual body RIRF matrix
    return file_info.GetRIRFVal(b, r, col, st);
}

Eigen::VectorXd TestHydro::ComputeForceWaves() {
    Eigen::VectorXd wf(num_bodies * 6);
    // TODO: check size of force calculated against wf (this could help catch/fix NoWave multibody issue)
    wf          = user_waves->GetForceAtTime(bodies[0]->GetChTime());
    force_waves = wf;
    return force_waves;
}

double TestHydro::coordinateFunc(int b, int i) {
    int body_num_offset = 6 * (b - 1);  // b_num from ForceFunc6d is 1 indexed, TODO: make all b_num 0 indexed
    int total_dofs      = 6 * num_bodies;
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
    // TODO change forces to Eigen::VectorXd types, might need to force evaluation of eigen at some point?
    std::fill(total_force.begin(), total_force.end(), 0.0);
    std::fill(force_hydrostatic.begin(), force_hydrostatic.end(), 0.0);
    std::fill(force_radiation_damping.begin(), force_radiation_damping.end(), 0.0);
    std::fill(force_waves.begin(), force_waves.end(), 0.0);

    // call compute forces
    convTrapz = true;  // use trapeziodal rule or assume fixed dt.

    force_hydrostatic       = ComputeForceHydrostatics();
    force_radiation_damping = ComputeForceRadiationDampingConv();  // TODO non convolution option
    force_waves             = ComputeForceWaves();

    // TODO once all force components are Eigen, remove this from being a loop and add vectors directly
    for (int i = 0; i < total_dofs; i++) {
        total_force[i] = force_hydrostatic[i] - force_radiation_damping[i] + force_waves[i];
    }

    // std::cout << "force_waves\n";
    // for (int i = 0; i < total_dofs; i++) {
    //    std::cout << force_waves[i] << std::endl;
    //}

    if (body_num_offset + i < 0 || body_num_offset >= total_dofs) {
        std::cout << "total force accessing out of bounds" << std::endl;
    }
    return total_force[body_num_offset + i];
}
