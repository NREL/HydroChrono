/*********************************************************************
 * @file  hydrostatics_model.cpp
 * @brief Implementation of hydrostatics force model.
 *********************************************************************/

#include "hydrostatics_model.h"
#include <hydroc/h5fileinfo.h>

#include <Eigen/Dense>
#include <cassert>

namespace hydrochrono::hydro {

HydrostaticsModel::HydrostaticsModel(
    const HydroData& file_info,
    int num_bodies,
    const std::vector<double>& equilibrium,
    const std::vector<double>& cb_minus_cg,
    const Eigen::Vector3d& gravitational_acceleration)
    : file_info_(file_info),
      num_bodies_(num_bodies),
      equilibrium_(equilibrium),
      cb_minus_cg_(cb_minus_cg),
      gravitational_acceleration_(gravitational_acceleration) {
    assert(num_bodies_ > 0);
}

void HydrostaticsModel::Compute(const SystemState& state,
                                double time,
                                BodyForces& inout_forces) {
    assert(static_cast<int>(state.bodies.size()) == num_bodies_);
    assert(static_cast<int>(inout_forces.size()) == num_bodies_);

    const double rho = file_info_.GetRhoVal();
    const double rho_times_g = rho * gravitational_acceleration_.norm();

    // Ensure forces are properly sized (6 DOF per body)
    for (int b = 0; b < num_bodies_; ++b) {
        if (static_cast<int>(inout_forces[b].size()) != kDofPerBody) {
            inout_forces[b].resize(kDofPerBody);
            inout_forces[b].setZero();
        }
    }

    // Multibody loop: compute hydrostatic forces for each body
    for (int b = 0; b < num_bodies_; ++b) {
        const auto& body_state = state.bodies[b];
        const int body_offset = kDofPerBody * b;
        const double* const body_equilibrium = &equilibrium_[body_offset];

        // Current pose (from SystemState)
        const Eigen::Vector3d& position_world = body_state.position;
        const Eigen::Vector3d& rotation_rpy = body_state.orientation_rpy;

        // 6-DOF displacement from equilibrium (translation xyz, rotation rpy)
        Eigen::Matrix<double, kDofPerBody, 1> displacement_from_equilibrium;
        displacement_from_equilibrium[0] = position_world[0] - body_equilibrium[0];
        displacement_from_equilibrium[1] = position_world[1] - body_equilibrium[1];
        displacement_from_equilibrium[2] = position_world[2] - body_equilibrium[2];
        displacement_from_equilibrium[3] = rotation_rpy[0] - body_equilibrium[3];
        displacement_from_equilibrium[4] = rotation_rpy[1] - body_equilibrium[4];
        displacement_from_equilibrium[5] = rotation_rpy[2] - body_equilibrium[5];

        // Linear hydrostatic restoring force/torque
        const Eigen::MatrixXd restoring_stiffness_matrix = file_info_.GetLinMatrix(b);
        const Eigen::Matrix<double, kDofPerBody, 1> restoring_force_torque =
            -rho_times_g * (restoring_stiffness_matrix * displacement_from_equilibrium);
        for (int i = 0; i < kDofPerBody; ++i) {
            inout_forces[b][i] += restoring_force_torque[i];
        }

        // Buoyancy force at equilibrium: F = rho * (-g) * displaced_volume
        const double displaced_volume = file_info_.GetDispVolVal(b);
        const Eigen::Vector3d buoyancy_force = rho * (-gravitational_acceleration_) * displaced_volume;
        inout_forces[b][0] += buoyancy_force[0];
        inout_forces[b][1] += buoyancy_force[1];
        inout_forces[b][2] += buoyancy_force[2];

        // Buoyancy-induced moment about CG: (r_CB - r_CG) x buoyancy
        const int rotation_offset = kDofLinOrRot * b;
        const Eigen::Vector3d cg_to_cb(
            cb_minus_cg_[rotation_offset + 0],
            cb_minus_cg_[rotation_offset + 1],
            cb_minus_cg_[rotation_offset + 2]
        );
        const Eigen::Vector3d buoyancy_torque = cg_to_cb.cross(buoyancy_force);
        inout_forces[b][3] += buoyancy_torque[0];
        inout_forces[b][4] += buoyancy_torque[1];
        inout_forces[b][5] += buoyancy_torque[2];
    }
}

}  // namespace hydrochrono::hydro

