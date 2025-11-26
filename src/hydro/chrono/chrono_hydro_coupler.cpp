/*********************************************************************
 * @file  chrono_hydro_coupler.cpp
 * @brief Implementation of ChronoHydroCoupler.
 *********************************************************************/

#include "chrono_hydro_coupler.h"

#include <cassert>

namespace hydrochrono::hydro {

ChronoHydroCoupler::ChronoHydroCoupler(
    std::shared_ptr<HydroSystem> hydro_system,
    std::vector<std::shared_ptr<chrono::ChBody>> bodies)
    : hydro_system_(hydro_system),
      bodies_(std::move(bodies)) {
    assert(hydro_system_ != nullptr);
    assert(!bodies_.empty());
}

BodyForces ChronoHydroCoupler::Evaluate(double time) const {
    // Build SystemState from Chrono bodies
    SystemState system_state;
    BuildSystemStateFromChronoBodies(bodies_, system_state);

    // Evaluate forces via HydroSystem
    return hydro_system_->Evaluate(system_state, time);
}

void ChronoHydroCoupler::ApplyForcesToChrono(const BodyForces& forces) {
    // TODO: Implement force application to Chrono bodies
    // This will push forces back into the Chrono physics engine
    (void)forces;  // Suppress unused parameter warning
}

}  // namespace hydrochrono::hydro

