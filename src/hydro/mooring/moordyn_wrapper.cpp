/*********************************************************************
 * @file  moordyn_wrapper.cpp
 * @brief Implementation of the MoorDyn RAII wrapper.
 *********************************************************************/

#include "moordyn_wrapper.h"
#include <MoorDyn2.h>

#include <stdexcept>
#include <string>
#include <cmath>
#include <cassert>

namespace hydrochrono::hydro {

static constexpr int kDofPerBody = 6;

MoorDynWrapper::MoorDynWrapper(const std::string& input_file,
                               const std::vector<int>& coupled_body_indices)
    : input_file_(input_file),
      coupled_body_indices_(coupled_body_indices) {
    if (input_file_.empty()) {
        throw std::invalid_argument(
            "MoorDynWrapper: input_file must not be empty");
    }
    if (coupled_body_indices_.empty()) {
        throw std::invalid_argument(
            "MoorDynWrapper: at least one coupled body index is required");
    }
}

MoorDynWrapper::~MoorDynWrapper() {
    if (system_) {
        MoorDyn_Close(system_);
        system_ = nullptr;
    }
}

void MoorDynWrapper::Initialize(const SystemState& initial_state) {
    if (initialized_) {
        throw std::runtime_error(
            "MoorDynWrapper::Initialize called more than once");
    }

    system_ = MoorDyn_Create(input_file_.c_str());
    if (!system_) {
        throw std::runtime_error(
            "MoorDyn_Create failed for input file: " + input_file_);
    }

    int err = MoorDyn_NCoupledDOF(system_, &n_coupled_dof_);
    if (err != 0) {
        throw std::runtime_error("MoorDyn_NCoupledDOF failed");
    }

    const unsigned int expected_dof =
        static_cast<unsigned int>(coupled_body_indices_.size()) * kDofPerBody;
    if (n_coupled_dof_ != expected_dof) {
        throw std::runtime_error(
            "MoorDyn reports " + std::to_string(n_coupled_dof_) +
            " coupled DOFs but HydroChrono expected " +
            std::to_string(expected_dof) + " (" +
            std::to_string(coupled_body_indices_.size()) +
            " bodies x 6 DOF). Check the MoorDyn input file body list.");
    }

    x_.resize(n_coupled_dof_, 0.0);
    xd_.resize(n_coupled_dof_, 0.0);
    f_.resize(n_coupled_dof_, 0.0);

    PackState(initial_state, x_, xd_);

    err = MoorDyn_Init(system_, x_.data(), xd_.data());
    if (err != 0) {
        throw std::runtime_error("MoorDyn_Init failed with error code " +
                                 std::to_string(err));
    }

    moordyn_time_ = 0.0;
    initialized_ = true;
}

void MoorDynWrapper::Step(const SystemState& state, double time, double dt,
                          BodyForces& out_forces) {
    assert(initialized_ && "MoorDynWrapper::Step called before Initialize");

    if (dt <= 0.0) {
        return;
    }

    PackState(state, x_, xd_);

    double md_dt = dt;
    int err = MoorDyn_Step(system_, x_.data(), xd_.data(), f_.data(),
                           &moordyn_time_, &md_dt);
    if (err != 0) {
        throw std::runtime_error(
            "MoorDyn_Step failed at t=" + std::to_string(time) +
            " with error code " + std::to_string(err));
    }

    for (size_t ci = 0; ci < coupled_body_indices_.size(); ++ci) {
        const int bi = coupled_body_indices_[ci];
        assert(bi >= 0 && bi < static_cast<int>(out_forces.size()));

        if (out_forces[bi].size() < kDofPerBody) {
            out_forces[bi].setZero(kDofPerBody);
        }

        const size_t offset = ci * kDofPerBody;
        for (int d = 0; d < kDofPerBody; ++d) {
            out_forces[bi](d) += f_[offset + d];
        }
    }
}

void MoorDynWrapper::PackState(const SystemState& state,
                               std::vector<double>& x,
                               std::vector<double>& xd) const {
    for (size_t ci = 0; ci < coupled_body_indices_.size(); ++ci) {
        const int bi = coupled_body_indices_[ci];
        assert(bi >= 0 && bi < static_cast<int>(state.bodies.size()));
        const BodyState& bs = state.bodies[bi];

        const size_t offset = ci * kDofPerBody;

        // Position: [x, y, z, roll, pitch, yaw]
        x[offset + 0] = bs.position[0];
        x[offset + 1] = bs.position[1];
        x[offset + 2] = bs.position[2];
        x[offset + 3] = bs.orientation_rpy[0];
        x[offset + 4] = bs.orientation_rpy[1];
        x[offset + 5] = bs.orientation_rpy[2];

        // Velocity: [vx, vy, vz, wx, wy, wz]
        xd[offset + 0] = bs.linear_velocity[0];
        xd[offset + 1] = bs.linear_velocity[1];
        xd[offset + 2] = bs.linear_velocity[2];
        xd[offset + 3] = bs.angular_velocity[0];
        xd[offset + 4] = bs.angular_velocity[1];
        xd[offset + 5] = bs.angular_velocity[2];
    }
}

}  // namespace hydrochrono::hydro
