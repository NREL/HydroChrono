/*********************************************************************
 * @file  hydro_system.cpp
 * @brief Implementation of HydroSystem façade.
 *********************************************************************/

#include <hydroc/core/hydro_system.h>

#include <chrono>
#include <stdexcept>
#include <Eigen/Dense>

namespace hydrochrono::hydro {

HydroSystem::HydroSystem(int num_bodies,
                         std::vector<std::unique_ptr<IHydroForceComponent>> components)
    : num_bodies_(num_bodies),
      components_(std::move(components)),
      profile_stats_(),
      profiling_enabled_(false) {
    // Require at least one body; otherwise hydrodynamic forces are meaningless.
    if (num_bodies_ <= 0) {
        throw std::invalid_argument(
            "HydroSystem: num_bodies must be > 0 (got " + std::to_string(num_bodies_) + ")");
    }
}

BodyForces HydroSystem::Evaluate(const SystemState& state, double time) {
    // Initialize forces: one 6-DOF zero vector per body
    BodyForces forces(num_bodies_);
    for (int b = 0; b < num_bodies_; ++b) {
        forces[b].resize(6);
        forces[b].setZero();
    }

    // Accumulate force contributions from all components
    for (const auto& component : components_) {
        if (profiling_enabled_) {
            // Profiling path: measure timing for each component
            auto t0 = std::chrono::steady_clock::now();
            component->Compute(state, time, forces);
            auto t1 = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(t1 - t0).count();
            
            switch (component->Type()) {
                case HydroComponentType::Hydrostatics:
                    profile_stats_.hydrostatics_seconds += elapsed;
                    profile_stats_.hydrostatics_calls++;
                    break;
                case HydroComponentType::Radiation:
                    profile_stats_.radiation_seconds += elapsed;
                    profile_stats_.radiation_calls++;
                    break;
                case HydroComponentType::Excitation:
                    profile_stats_.excitation_seconds += elapsed;
                    profile_stats_.excitation_calls++;
                    break;
            }
        } else {
            // Fast path: no timing overhead
            component->Compute(state, time, forces);
        }
    }

    return forces;
}

}  // namespace hydrochrono::hydro

