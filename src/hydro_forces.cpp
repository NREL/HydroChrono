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

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>  // std::accumulate
#include <random>
#include <stdexcept>
#include <vector>

const int kDofPerBody  = 6;
const int kDofLinOrRot = 3;

/**
 * @brief Generates a vector of evenly spaced numbers over a specified range.
 *
 * This function returns a vector of `num_points` numbers evenly spaced from
 * `start` to `end`. The function utilizes a single loop for this computation,
 * making it efficient for generating large vectors.
 *
 * @param start - The start value of the sequence.
 * @param end - The end value of the sequence.
 * @param num_points - The number of evenly spaced samples to generate.
 * @return std::vector<double> - Vector of evenly spaced numbers.
 * @exception None
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
    index_ = kDofPerBody;
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
    std::string temp  = body_->GetNameString();  // remove "body" from "bodyN", convert N to int, get body num
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
    if (i >= kDofPerBody || i < 0) {
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

TestHydro::TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
                     std::string h5_file_name,
                     std::shared_ptr<WaveBase> waves)
    : bodies_(user_bodies),
      num_bodies_(bodies_.size()),
      file_info_(H5FileInfo(h5_file_name, num_bodies_).ReadH5Data()) {
    prev_time = -1;

    // Set up time vector
    rirf_time_vector = file_info_.GetRIRFTimeVector();
    // width array
    rirf_width_vector.resize(rirf_time_vector.size());
    for (int ii = 0; ii < rirf_width_vector.size(); ii++) {
        rirf_width_vector[ii] = 0.0;
        if (ii < rirf_time_vector.size() - 1) {
            rirf_width_vector[ii] += 0.5 * abs(rirf_time_vector[ii + 1] - rirf_time_vector[ii]);
        }
        if (ii > 0) {
            rirf_width_vector[ii] += 0.5 * abs(rirf_time_vector[ii] - rirf_time_vector[ii - 1]);
        }
    }

    // Total degrees of freedom
    int total_dofs = kDofPerBody * num_bodies_;

    // Initialize vectors
    time_history_.clear();
    velocity_history_.clear();
    for (int b = 0; b < num_bodies_; ++b) {
        velocity_history_.push_back(std::vector<std::vector<double>>(0));
    }
    force_hydrostatic_.assign(total_dofs, 0.0);
    force_radiation_damping_.assign(total_dofs, 0.0);
    total_force_.assign(total_dofs, 0.0);
    equilibrium_.assign(total_dofs, 0.0);
    cb_minus_cg_.assign(kDofLinOrRot * num_bodies_, 0.0);

    // Compute equilibrium and cb_minus_cg_
    for (int b = 0; b < num_bodies_; ++b) {
        for (int i = 0; i < kDofLinOrRot; ++i) {
            unsigned eq_idx = i + kDofPerBody * b;
            unsigned c_idx  = i + kDofLinOrRot * b;

            equilibrium_[eq_idx] = file_info_.GetCGVector(b)[i];
            cb_minus_cg_[c_idx]  = file_info_.GetCBVector(b)[i] - file_info_.GetCGVector(b)[i];
        }
    }

    for (int b = 0; b < num_bodies_; ++b) {
        force_per_body_.emplace_back(bodies_[b], this);
    }

    // Handle added mass info
    my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();

    std::vector<std::shared_ptr<ChLoadable>> loadables(bodies_.size());
    for (int i = 0; i < static_cast<int>(bodies_.size()); ++i) {
        loadables[i] = bodies_[i];
    }

    my_loadbodyinertia =
        chrono_types::make_shared<ChLoadAddedMass>(file_info_.GetBodyInfos(), loadables, bodies_[0]->GetSystem());

    bodies_[0]->GetSystem()->Add(my_loadcontainer);
    my_loadcontainer->Add(my_loadbodyinertia);

    // Set up hydro inputs
    user_waves_ = waves;
    AddWaves(user_waves_);
}

void TestHydro::AddWaves(std::shared_ptr<WaveBase> waves) {
    user_waves_ = waves;

    switch (user_waves_->GetWaveMode()) {
        case WaveMode::regular: {
            auto reg = std::static_pointer_cast<RegularWave>(user_waves_);
            reg->AddH5Data(file_info_.GetRegularWaveInfos(), file_info_.GetSimulationInfo());
            break;
        }
        case WaveMode::irregular: {
            auto irreg = std::static_pointer_cast<IrregularWaves>(user_waves_);
            irreg->AddH5Data(file_info_.GetIrregularWaveInfos(), file_info_.GetSimulationInfo());
            break;
        }
    }

    user_waves_->Initialize();
}

std::vector<double> TestHydro::ComputeForceHydrostatics() {
    assert(num_bodies_ > 0);

    const double rho = file_info_.GetRhoVal();
    const auto g_acc = bodies_[0]->GetSystem()->Get_G_acc();  // assuming all bodies in same system
    const double gg  = g_acc.Length();

    for (int b = 0; b < num_bodies_; b++) {
        std::shared_ptr<chrono::ChBody> body = bodies_[b];

        int b_offset                   = kDofPerBody * b;
        double* body_force_hydrostatic = &force_hydrostatic_[b_offset];
        double* body_equilibrium       = &equilibrium_[b_offset];

        // hydrostatic stiffness due to offset from equilibrium
        const auto body_position = body->GetPos();
        const auto body_rotation = body->GetRot().Q_to_Euler123();

        chrono::ChVectorN<double, kDofPerBody> body_displacement;
        for (int ii = 0; ii < kDofLinOrRot; ii++) {
            body_displacement[ii]                = body_position[ii] - body_equilibrium[ii];
            body_displacement[ii + kDofLinOrRot] = body_rotation[ii] - body_equilibrium[ii + kDofLinOrRot];
        }

        const auto force_offset = -gg * rho * file_info_.GetLinMatrix(b) * body_displacement;
        for (int dof = 0; dof < kDofPerBody; dof++) {
            body_force_hydrostatic[dof] += force_offset[dof];
        }

        // buoyancy at equilibrium
        const auto buoyancy = rho * (-g_acc) * file_info_.GetDispVolVal(b);

        for (int ii = 0; ii < kDofLinOrRot; ii++) {
            body_force_hydrostatic[ii] += buoyancy[ii];
        }

        int r_offset = kDofLinOrRot * b;
        const auto cg2cb =
            chrono::ChVector<double>(cb_minus_cg_[r_offset], cb_minus_cg_[r_offset + 1], cb_minus_cg_[r_offset + 2]);
        const auto buoyancy2 = cg2cb % buoyancy;

        for (int ii = 0; ii < kDofLinOrRot; ii++) {
            body_force_hydrostatic[ii + kDofLinOrRot] += buoyancy2[ii];
        }
    }

    return force_hydrostatic_;
}

std::vector<double> TestHydro::ComputeForceRadiationDampingConv() {
    const int size    = file_info_.GetRIRFDims(2);
    const int numRows = kDofPerBody * num_bodies_;
    const int numCols = kDofPerBody * num_bodies_;

    assert(numRows * size > 0 && numCols > 0);

    // time history
    auto t_sim = bodies_[0]->GetChTime();
    auto t_min = t_sim - rirf_time_vector.tail<1>()[0];
    if (time_history_.size() > 0 && t_sim == time_history_.front()) {
        throw std::runtime_error("Tried to compute the radiation damping convolution twice within the same time step!");
    }
    time_history_.insert(time_history_.begin(), t_sim);

    // velocity history
    for (int b = 0; b < num_bodies_; b++) {
        auto& body                  = bodies_[b];
        auto& velocity_history_body = velocity_history_[b];

        auto vel                    = body->GetPos_dt();
        auto wvel                   = body->GetWvel_par();
        std::vector<double> vel_vec = {vel[0], vel[1], vel[2], wvel[0], wvel[1], wvel[2]};
        velocity_history_body.insert(velocity_history_body.begin(), vel_vec);
    }

    // remove unnecessary history
    if (time_history_.size() > 1) {
        while (time_history_[time_history_.size() - 2] < t_min) {
            time_history_.pop_back();
            for (int b = 0; b < num_bodies_; b++) {
                auto& velocity_history_body = velocity_history_[b];
                velocity_history_body.pop_back();
            }
        }
    }

    if (time_history_.size() > 1) {
        int idx_history = 0;

        // iterate over RIRF steps
        for (int step = 0; step < size; step++) {
            auto t_rirf = t_sim - rirf_time_vector[step];
            while (time_history_[idx_history + 1] > t_rirf && idx_history < time_history_.size() - 1) {
                idx_history += 1;
            }
            if (idx_history >= time_history_.size() - 1) {
                break;
            }

            // iterate over bodies
            for (int idx_body = 0; idx_body < num_bodies_; idx_body++) {
                auto& velocity_history_body = velocity_history_[idx_body];
                std::vector<double> vel(kDofPerBody);

                // interpolate velocity at t_rirf from recorded velocity history
                // time values
                auto t1 = time_history_[idx_history + 1];
                auto t2 = time_history_[idx_history];
                if (t_rirf == t1) {
                    vel = velocity_history_body[idx_history + 1];
                } else if (t_rirf == t2) {
                    vel = velocity_history_body[idx_history];
                } else if (t_rirf > t1 && t_rirf < t2) {
                    // weights
                    auto w1 = (t2 - t_rirf) / (t2 - t1);
                    auto w2 = 1.0 - w1;
                    // velocity values
                    auto vel1 = velocity_history_body[idx_history + 1];
                    auto vel2 = velocity_history_body[idx_history];
                    for (int dof = 0; dof < kDofPerBody; dof++) {
                        vel[dof] = w1 * vel1[dof] + w2 * vel2[dof];
                    }
                } else {
                    throw std::runtime_error("Radiation convolution: wrong interpolation: " + std::to_string(t_rirf) +
                                             " not between " + std::to_string(t1) + " and " + std::to_string(t2) + ".");
                }

                for (int dof = 0; dof < kDofPerBody; dof++) {
                    // get column index
                    int col = dof + idx_body * kDofPerBody;

                    // iterate over rows
                    for (int row = 0; row < numRows; row++) {
                        force_radiation_damping_[row] +=
                            GetRIRFval(row, col, step) * vel[dof] * rirf_width_vector[step];
                    }
                }
            }
        }
    }
    return force_radiation_damping_;
}

double TestHydro::GetRIRFval(int row, int col, int st) {
    if (row < 0 || row >= kDofPerBody * num_bodies_ || col < 0 || col >= kDofPerBody * num_bodies_ || st < 0 ||
        st >= file_info_.GetRIRFDims(2)) {
        throw std::out_of_range("rirfval index out of range in TestHydro");
    }

    int body_index = row / kDofPerBody;
    int col_dof    = col % kDofPerBody;
    int row_dof    = row % kDofPerBody;

    return file_info_.GetRIRFVal(body_index, row_dof, col, st);
}

Eigen::VectorXd TestHydro::ComputeForceWaves() {
    // Ensure bodies_ is not empty
    if (bodies_.empty()) {
        throw std::runtime_error("bodies_ array is empty in ComputeForceWaves");
    }

    force_waves_ = user_waves_->GetForceAtTime(bodies_[0]->GetChTime());

    // TODO: Add size check for force_waves_ if needed
    // Example:
    // if (force_waves_.size() != expected_size) {
    //     throw std::runtime_error("Mismatched size in ComputeForceWaves");
    // }

    return force_waves_;
}

double TestHydro::CoordinateFuncForBody(int b, int dof_index) {
    if (dof_index < 0 || dof_index >= kDofPerBody || b < 1 || b > num_bodies_) {
        throw std::out_of_range("Invalid index in CoordinateFuncForBody");
    }

    // Adjusting for 1-indexed body number
    const int body_num_offset = kDofPerBody * (b - 1);
    const int total_dofs      = kDofPerBody * num_bodies_;

    // Ensure the bodies_ vector isn't empty and the first element isn't null
    if (bodies_.empty() || !bodies_[0]) {
        throw std::runtime_error("bodies_ array is empty or invalid in CoordinateFuncForBody");
    }

    // Check if the forces for this time step have already been computed
    if (bodies_[0]->GetChTime() == prev_time) {
        return total_force_[body_num_offset + dof_index];
    }

    // Update time and reset forces for this time step
    prev_time = bodies_[0]->GetChTime();
    std::fill(total_force_.begin(), total_force_.end(), 0.0);
    std::fill(force_hydrostatic_.begin(), force_hydrostatic_.end(), 0.0);
    std::fill(force_radiation_damping_.begin(), force_radiation_damping_.end(), 0.0);
    std::fill(force_waves_.begin(), force_waves_.end(), 0.0);

    force_hydrostatic_       = ComputeForceHydrostatics();
    force_radiation_damping_ = ComputeForceRadiationDampingConv();
    force_waves_             = ComputeForceWaves();

    // Accumulate total force (consider converting forces to Eigen::VectorXd in the future for direct addition)
    for (int index = 0; index < total_dofs; index++) {
        total_force_[index] = force_hydrostatic_[index] - force_radiation_damping_[index] + force_waves_[index];
    }

    if (body_num_offset + dof_index < 0 || body_num_offset >= total_dofs) {
        throw std::out_of_range("Accessing out-of-bounds index in CoordinateFuncForBody");
    }

    return total_force_[body_num_offset + dof_index];
}