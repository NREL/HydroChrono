/*********************************************************************
 * @file  excitation_model.cpp
 * @brief Implementation of wave excitation force model.
 *********************************************************************/

#include "excitation_model.h"

#include <cassert>
#include <Eigen/Dense>

namespace hydrochrono::hydro {

ExcitationModel::ExcitationModel(std::shared_ptr<WaveBase> waves, int num_bodies)
    : waves_(waves), num_bodies_(num_bodies) {
    assert(waves_ != nullptr);
    assert(num_bodies_ > 0);
}

void ExcitationModel::Compute(const SystemState& state,
                              double time,
                              BodyForces& inout_forces) {
    assert(static_cast<int>(inout_forces.size()) == num_bodies_);

    // Get wave excitation forces from WaveBase (returns flat 6N vector)
    Eigen::VectorXd force_flat = waves_->GetForceAtTime(time);

    // Ensure forces are properly sized (6 DOF per body)
    for (int b = 0; b < num_bodies_; ++b) {
        if (static_cast<int>(inout_forces[b].size()) != kDofPerBody) {
            inout_forces[b].resize(kDofPerBody);
            inout_forces[b].setZero();
        }
    }

    // Convert flat 6N vector to per-body forces and add to inout_forces
    for (int b = 0; b < num_bodies_; ++b) {
        const int body_offset = kDofPerBody * b;
        for (int i = 0; i < kDofPerBody; ++i) {
            inout_forces[b][i] += force_flat[body_offset + i];
        }
    }
}

}  // namespace hydrochrono::hydro

