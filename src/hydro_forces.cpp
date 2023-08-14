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
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <vector>


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
    base_  = NULL;
    index_ = 6;
}

ComponentFunc::ComponentFunc(ForceFunc6d* b, int i) : base_(b), index_(i) {}

ComponentFunc* ComponentFunc::Clone() const {
    return new ComponentFunc(*this);
}

ComponentFunc::ComponentFunc(const ComponentFunc& old) {
    base_  = old.base_;
    index_ = old.index_;
}

double ComponentFunc::Get_y(double x) const {
    if (base_ == NULL) {
        std::cout << "base == Null!" << std::endl;
        return 0;
    }
    return base_->CoordinateFunc(index_);
}

ForceFunc6d::ForceFunc6d() : forces_{{this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5}} {
    for (unsigned i = 0; i < 6; i++) {
        force_ptrs_[i] = std::shared_ptr<ComponentFunc>(forces_ + i, [](ComponentFunc*) {});
        // sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is
        // automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
        // them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the
        // default deletion logic to do nothing
        // Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
    }
    chrono_force_  = chrono_types::make_shared<ChForce>();
    chrono_torque_ = chrono_types::make_shared<ChForce>();
    chrono_force_->SetAlign(ChForce::AlignmentFrame::WORLD_DIR);
    chrono_torque_->SetAlign(ChForce::AlignmentFrame::WORLD_DIR);
    chrono_force_->SetNameString("hydroforce");
    chrono_torque_->SetNameString("hydrotorque");
}

ForceFunc6d::ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* user_all_forces) : ForceFunc6d() {
    body_             = object;
    std::string temp = body_->GetNameString();   // remove "body" from "bodyN", convert N to int, get body num
    b_num_            = stoi(temp.erase(0, 4));  // 1 indexed TODO: fix b_num starting here to be 0 indexed
    all_hydro_forces_ = user_all_forces;         // TODO switch to smart pointers? does this use = ?
    if (all_hydro_forces_ == NULL) {
        std::cout << "all hydro forces null " << std::endl;
    }
    SetForce();
    SetTorque();
    ApplyForceAndTorqueToBody();
}

ForceFunc6d::ForceFunc6d(const ForceFunc6d& old)
    : forces_{{this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5}} {
    for (unsigned i = 0; i < 6; i++) {
        force_ptrs_[i] = std::shared_ptr<ComponentFunc>(forces_ + i, [](ComponentFunc*) {});
        // sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is
        // automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
        // them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the
        // default deletion logic to do nothing
        // Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
    }
    chrono_force_     = old.chrono_force_;
    chrono_torque_    = old.chrono_torque_;
    body_             = old.body_;
    b_num_            = old.b_num_;
    all_hydro_forces_ = old.all_hydro_forces_;
    SetForce();
    SetTorque();
}

double ForceFunc6d::CoordinateFunc(int i) {
    // b_num is 1 indexed?
    if (i >= 6 || i < 0) {
        std::cout << "wrong index force func 6d" << std::endl;
        return 0;
    }
    return all_hydro_forces_->CoordinateFuncForBody(
        b_num_, i);  // b_num is 1 indexed here!!!!! TODO: change all b_num to be 0 indexed everywhere
}

void ForceFunc6d::SetForce() {
    if (chrono_force_ == NULL || body_ == NULL) {
        std::cout << "set force null issue" << std::endl;
    }
    chrono_force_->SetF_x(force_ptrs_[0]);
    chrono_force_->SetF_y(force_ptrs_[1]);
    chrono_force_->SetF_z(force_ptrs_[2]);
}

void ForceFunc6d::SetTorque() {
    if (chrono_torque_ == NULL || body_ == NULL) {
        std::cout << "set torque null issue" << std::endl;
    }
    chrono_torque_->SetF_x(force_ptrs_[3]);
    chrono_torque_->SetF_y(force_ptrs_[4]);
    chrono_torque_->SetF_z(force_ptrs_[5]);
    chrono_torque_->SetMode(ChForce::ForceType::TORQUE);
}

void ForceFunc6d::ApplyForceAndTorqueToBody() {
    body_->AddForce(chrono_force_);
    body_->AddForce(chrono_torque_);
}

// TODO reorder TestHydro function ordering to match header file order
TestHydro::TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
                     std::string h5_file_name,
                     std::shared_ptr<WaveBase> waves)
    : bodies_(user_bodies), num_bodies_(bodies_.size()), file_info_(H5FileInfo(h5_file_name, num_bodies_).ReadH5Data()) {
    prev_time   = -1;
    offset_rirf = 0;

    // set up time vector (should be the same for each body, so just use the first always)
    rirf_time_vector = file_info_.GetRIRFTimeVector();
    rirf_timestep_    = rirf_time_vector[1] - rirf_time_vector[0];  // TODO is this the same for all bodies?

    // simplify 6* num_bodies to be the system's total number of dofs, makes expressions later easier to read
    int total_dofs = 6 * num_bodies_;
    // resize and initialize velocity history vector to all zeros
    velocity_history.resize(file_info_.GetRIRFDims(2) * total_dofs, 0.0);  // resize and fill with 0s
    // resize and initialize all persistent forces to all 0s
    // TODO rephrase for Eigen::VectorXd eventually
    force_hydrostatic_.resize(total_dofs, 0.0);
    force_radiation_damping_.resize(total_dofs, 0.0);
    total_force_.resize(total_dofs, 0.0);
    // set up equilibrium for entire system (each body has position and rotation equilibria 3 indicies apart)
    equilibrium_.resize(total_dofs, 0.0);
    cb_minus_cg_.resize(3 * num_bodies_, 0.0);  // cb-cg has 3 components for each body
    for (int b = 0; b < num_bodies_; b++) {
        for (int i = 0; i < 3; i++) {
            unsigned equilibrium_idx = i + 6 * b;
            unsigned c_idx           = i + 3 * b;
            // positional equilib is cg, leave rotational bit 0
            equilibrium_[equilibrium_idx] = file_info_.GetCGVector(b)[i];
            cb_minus_cg_[c_idx]           = file_info_.GetCBVector(b)[i] - file_info_.GetCGVector(b)[i];
        }
    }

    for (int b = 0; b < num_bodies_; b++) {
        force_per_body_.emplace_back(bodies_[b], this);
    }

    // added mass info
    my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();

    /// TODO Check if local vector is really copied into constructor of ChLoadAddedMass
    /// else it could be a memory fault
    std::vector<std::shared_ptr<ChLoadable>> loadables(bodies_.size());
    for (auto i = 0; i < bodies_.size(); ++i) {
        loadables[i] = bodies_[i];
    }

    my_loadbodyinertia =
        chrono_types::make_shared<ChLoadAddedMass>(file_info_.GetBodyInfos(), loadables, bodies_[0]->GetSystem());
    bodies_[0]->GetSystem()->Add(my_loadcontainer);
    my_loadcontainer->Add(my_loadbodyinertia);

    // set up hydro inputs stuff
    // hydro_inputs = user_hydro_inputs;
    // WaveSetUp();
    user_waves_ = waves;
    AddWaves(user_waves_);
}

void TestHydro::AddWaves(std::shared_ptr<WaveBase> waves) {
    user_waves_ = waves;
    if (user_waves_->GetWaveMode() == WaveMode::regular) {
        std::shared_ptr<RegularWave> reg = std::static_pointer_cast<RegularWave>(user_waves_);
        reg->AddH5Data(file_info_.GetRegularWaveInfos());
    } else if (user_waves_->GetWaveMode() == WaveMode::irregular) {
        std::shared_ptr<IrregularWave> irreg = std::static_pointer_cast<IrregularWave>(user_waves_);
        irreg->AddH5Data(file_info_.GetIrregularWaveInfos(), file_info_.GetSimulationInfo());
    }
    user_waves_->Initialize();
}

double TestHydro::GetVelHistoryVal(int step, int c) const {
    if (step < 0 || step >= file_info_.GetRIRFDims(2) || c < 0 || c >= num_bodies_ * 6) {
        std::cout << "wrong vel history index " << std::endl;
        return 0;
    }
    int index = c % 6;
    int b     = c / 6;  // 0 indexed
    if (index + 6 * b + 6 * num_bodies_ * step >= num_bodies_ * 6 * file_info_.GetRIRFDims(2) ||
        index + 6 * b + 6 * num_bodies_ * step < 0) {
        std::cout << "bad vel history math" << std::endl;
        return 0;
    }
    return velocity_history[index + (6 * b) + (6 * num_bodies_ * step)];
}

double TestHydro::SetVelHistory(double val, int step, int b_num, int index) {
    if (step < 0 || step >= file_info_.GetRIRFDims(2) || b_num < 1 || b_num > num_bodies_ || index < 0 || index >= 6) {
        std::cout << "bad set vel history indexing" << std::endl;
        return 0;
    }
    if (index + 6 * (b_num - 1) + 6 * num_bodies_ * step < 0 ||
        index + 6 * (b_num - 1) + 6 * num_bodies_ * step >= num_bodies_ * 6 * file_info_.GetRIRFDims(2)) {
        std::cout << "bad set vel history math" << std::endl;
        return 0;
    }
    velocity_history[index + (6 * (b_num - 1)) + (6 * num_bodies_ * step)] = val;
    return val;
}

std::vector<double> TestHydro::ComputeForceHydrostatics() {
    assert(num_bodies_ > 0);

    for (int b = 0; b < num_bodies_; b++) {
        // initialize variables
        std::shared_ptr<chrono::ChBody> body = bodies_[b];
        // H5FileInfo& body_h5file              = file_info[b];
        double rho   = file_info_.GetRhoVal();
        int b_offset = 6 * b;
        // force_hydrostatic has 6 elements for each body so to skip to the next body we move 6 spaces
        double* body_force_hydrostatic = &force_hydrostatic_[b_offset];
        double* body_equilibrium       = &equilibrium_[b_offset];
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
        chrono::ChVectorN<double, 6> force_offset = -gg * rho * file_info_.GetLinMatrix(b) * body_displacement;
        // add to force_hydrostatic
        for (int dof = 0; dof < 6; dof++) {
            body_force_hydrostatic[dof] += force_offset[dof];
        }

        // buoyancy at equilibrium
        // TODO: move to prestep (shouldn't be calculated at each time step)
        // translational
        chrono::ChVector<> buoyancy =
            rho * (-body->GetSystem()->Get_G_acc()) * file_info_.GetDispVolVal(b);  // buoyancy = rho*g*Vdisp
        body_force_hydrostatic[0] += buoyancy[0];
        body_force_hydrostatic[1] += buoyancy[1];
        body_force_hydrostatic[2] += buoyancy[2];
        // rotational
        int r_offset = 3 * b;
        // cb_minus_cg has 3 elements for each body so to skip to the next body we move 3 spaces
        auto cg2cb =
            chrono::ChVector<double>(cb_minus_cg_[0 + r_offset], cb_minus_cg_[1 + r_offset], cb_minus_cg_[2 + r_offset]);
        chrono::ChVector<> buoyancy2 = cg2cb % buoyancy;
        body_force_hydrostatic[3] += buoyancy2[0];
        body_force_hydrostatic[4] += buoyancy2[1];
        body_force_hydrostatic[5] += buoyancy2[2];
    }
    return force_hydrostatic_;
}

std::vector<double> TestHydro::ComputeForceRadiationDampingConv() {
    int size = file_info_.GetRIRFDims(2);
    int nDoF = 6;
    // "shift" everything left 1
    offset_rirf--;  // starts as 0 before timestep change
    // keep offset close to 0, avoids small chance of -overflow errors in long simulations
    if (offset_rirf < -1 * size) {
        offset_rirf += size;
    }
    int numRows = nDoF * num_bodies_;
    int numCols = nDoF * num_bodies_;
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
        for (int b = 1; b < num_bodies_ + 1; b++) {  // body index being 1 indexed here is right
            int vi = (((size + offset_rirf) % size) + size) % size;
            SetVelHistory(bodies_[b - 1]->GetPos_dt()[i], vi, b, i);
            SetVelHistory(bodies_[b - 1]->GetWvel_par()[i], vi, b, i + 3);
        }
    }
    int vi;
    //#pragma omp parallel for
    if (convTrapz_ == true) {  // TODO add public function to TestHydro for users to change the value of convTrapz from
                              // main simulation
        // convolution integral using trapezoidal rule
        for (int row = 0; row < numRows; row++) {  // row goes to 6N
            for (int st = 0; st < size; st++) {
                vi = (((st + offset_rirf) % size) + size) % size;  // vi takes care of circshift function from matLab
                TMP_S(row, st) = 0;
                for (int col = 0; col < numCols; col++) {  // numCols goes to 6N
                    // multiply rirf by velocity history for each step and row (0,...,6N), store product in TIMESERIES
                    TIMESERIES(row, col, st) = GetRIRFval(row, col, st) * GetVelHistoryVal(vi, col);
                    // for (int i = 0; i < numCols; i++) {
                    // velOut << getVelHistoryVal(vi, col) << std::endl;
                    //}
                    // TMP_S is the sum over col (sum the effects of all radiating dofs (LDOF) for each time and motion
                    // dof)
                    TMP_S(row, st) += TIMESERIES(row, col, st);
                }
                if (st > 0) {
                    // integrate TMP_S
                    force_radiation_damping_[row] +=
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

    return force_radiation_damping_;
}

double TestHydro::GetRIRFval(int row, int col, int st) {
    if (row < 0 || row >= 6 * num_bodies_ || col < 0 || col >= 6 * num_bodies_ || st < 0 ||
        st >= file_info_.GetRIRFDims(2)) {
        std::cout << "rirfval index bad from testhydro" << std::endl;
        return 0;
    }
    int b = row / 6;  // 0 indexed, which body to get matrix info from
    int c = col % 6;  // which dof across column, 0,..,11 for 2 bodies, 0,...,6N-1 for N
    int r = row % 6;  // which dof 0,..,5 in individual body RIRF matrix
    return file_info_.GetRIRFVal(b, r, col, st);
}

Eigen::VectorXd TestHydro::ComputeForceWaves() {
    Eigen::VectorXd wf(num_bodies_ * 6);
    // TODO: check size of force calculated against wf (this could help catch/fix NoWave multibody issue)
    wf          = user_waves_->GetForceAtTime(bodies_[0]->GetChTime());
    force_waves_ = wf;
    return wf;
}

double TestHydro::CoordinateFuncForBody(int b, int i) {
    int body_num_offset = 6 * (b - 1);  // b_num from ForceFunc6d is 1 indexed, TODO: make all b_num 0 indexed
    int total_dofs      = 6 * num_bodies_;
    if (i < 0 || i > 5 || b < 1 || b > num_bodies_) {
        std::cout << "wrong index somewhere\nsetting coordinateFunc to 0" << std::endl;
        return 0;
    }
    // check prev_time here and only here
    // if forces have been computed for this time already, return the computed total force
    if (bodies_[0] == NULL) {
        std::cout << "bodies empty" << std::endl;
        return 0;
    }
    if (bodies_[0]->GetChTime() == prev_time) {
        return total_force_[body_num_offset + i];
    }
    // update current time and total_force for this step
    prev_time = bodies_[0]->GetChTime();

    // reset forces to 0
    // TODO change forces to Eigen::VectorXd types, might need to force evaluation of eigen at some point?
    std::fill(total_force_.begin(), total_force_.end(), 0.0);
    std::fill(force_hydrostatic_.begin(), force_hydrostatic_.end(), 0.0);
    std::fill(force_radiation_damping_.begin(), force_radiation_damping_.end(), 0.0);
    std::fill(force_waves_.begin(), force_waves_.end(), 0.0);

    // call compute forces
    convTrapz_ = true;  // use trapeziodal rule or assume fixed dt.

    force_hydrostatic_       = ComputeForceHydrostatics();
    force_radiation_damping_ = ComputeForceRadiationDampingConv();  // TODO non convolution option
    force_waves_             = ComputeForceWaves();

    // TODO once all force components are Eigen, remove this from being a loop and add vectors directly
    for (int i = 0; i < total_dofs; i++) {
        total_force_[i] = force_hydrostatic_[i] - force_radiation_damping_[i] + force_waves_[i];
    }

    // std::cout << "force_waves\n";
    // for (int i = 0; i < total_dofs; i++) {
    //    std::cout << force_waves[i] << std::endl;
    //}

    if (body_num_offset + i < 0 || body_num_offset >= total_dofs) {
        std::cout << "total force accessing out of bounds" << std::endl;
    }
    return total_force_[body_num_offset + i];
}
