/*********************************************************************
 * @file  excitation_component.cpp
 * @brief Implementation of wave excitation force component.
 *********************************************************************/

#include "excitation_component.h"

#include <stdexcept>
#include <Eigen/Dense>

namespace hydrochrono::hydro {

ExcitationComponent::ExcitationComponent(std::shared_ptr<WaveBase> waves, int num_bodies)
    : waves_(waves), num_bodies_(num_bodies) {
    // Wave model is required for excitation force computation.
    if (waves_ == nullptr) {
        throw std::invalid_argument(
            "ExcitationComponent: waves pointer must not be null");
    }

    // Require at least one body for wave excitation computation.
    if (num_bodies_ <= 0) {
        throw std::invalid_argument(
            "ExcitationComponent: num_bodies must be > 0 (got " + std::to_string(num_bodies_) + ")");
    }
}

void ExcitationComponent::Compute(const SystemState& state,
                              double time,
                              BodyForces& inout_forces) {
    // Internal consistency check: forces must match expected body count.
    if (static_cast<int>(inout_forces.size()) != num_bodies_) {
        throw std::runtime_error(
            "ExcitationComponent::Compute: inout_forces.size() mismatch (expected " +
            std::to_string(num_bodies_) + ", got " + std::to_string(inout_forces.size()) + ")");
    }

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

