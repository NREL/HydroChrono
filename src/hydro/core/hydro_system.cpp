/*********************************************************************
 * @file  hydro_system.cpp
 * @brief Implementation of HydroSystem façade.
 *********************************************************************/

#include "hydro_system.h"

#include <cassert>
#include <Eigen/Dense>

namespace hydrochrono::hydro {

HydroSystem::HydroSystem(int num_bodies,
                         std::vector<std::unique_ptr<IHydroForceComponent>> components)
    : num_bodies_(num_bodies),
      components_(std::move(components)) {
    assert(num_bodies_ > 0);
}

BodyForces HydroSystem::Evaluate(const SystemState& state, double time) const {
    // Initialize forces: one 6-DOF zero vector per body
    BodyForces forces(num_bodies_);
    for (int b = 0; b < num_bodies_; ++b) {
        forces[b].resize(6);
        forces[b].setZero();
    }

    // Accumulate force contributions from all components
    for (const auto& component : components_) {
        component->Compute(state, time, forces);
    }

    return forces;
}

}  // namespace hydrochrono::hydro

