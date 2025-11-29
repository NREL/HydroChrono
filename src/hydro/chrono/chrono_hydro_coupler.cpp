/*********************************************************************
 * @file  chrono_hydro_coupler.cpp
 * @brief Implementation of ChronoHydroCoupler.
 *********************************************************************/

#include <hydroc/coupling/chrono_coupler.h>
#include "chrono_state_utils.h"

#include <stdexcept>

namespace hydrochrono::hydro {

ChronoHydroCoupler::ChronoHydroCoupler(
    std::shared_ptr<HydroSystem> hydro_system,
    std::vector<std::shared_ptr<chrono::ChBody>> bodies)
    : hydro_system_(hydro_system),
      bodies_(std::move(bodies)) {
    // HydroSystem is required for force evaluation.
    if (hydro_system_ == nullptr) {
        throw std::invalid_argument(
            "ChronoHydroCoupler: hydro_system pointer must not be null");
    }

    // At least one body is required for hydrodynamic coupling.
    if (bodies_.empty()) {
        throw std::invalid_argument(
            "ChronoHydroCoupler: bodies vector must not be empty");
    }
}

BodyForces ChronoHydroCoupler::Evaluate(double time) {
    // Build SystemState from Chrono bodies
    SystemState system_state;
    chrono_coupling::BuildSystemStateFromChronoBodies(bodies_, system_state);

    // Evaluate forces via HydroSystem
    return hydro_system_->Evaluate(system_state, time);
}

void ChronoHydroCoupler::ApplyForcesToChrono(const BodyForces& forces) {
    // TODO: Implement force application to Chrono bodies
    // This will push forces back into the Chrono physics engine
    (void)forces;  // Suppress unused parameter warning
}

}  // namespace hydrochrono::hydro

