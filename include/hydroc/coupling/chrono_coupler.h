/*********************************************************************
 * @file  chrono_coupler.h
 * @brief ChronoHydroCoupler: bridges Chrono bodies and HydroForces.
 *
 * MAIN TYPES:
 *   - ChronoHydroCoupler: Adapter between Chrono and HydroForces
 *
 * ROLE: Extracts SystemState from Chrono bodies, invokes HydroForces
 * to compute forces. Used internally by HydroSystem façade.
 *********************************************************************/

#ifndef HYDROC_COUPLING_CHRONO_COUPLER_H
#define HYDROC_COUPLING_CHRONO_COUPLER_H

#include <hydroc/core/system_state.h>
#include <hydroc/core/hydro_forces.h>
#include <chrono/physics/ChBody.h>
#include <vector>
#include <memory>

namespace hydrochrono::hydro {

// Forward declaration for internal utilities
void BuildSystemStateFromChronoBodies(
    const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
    SystemState& out_state);

/**
 * @brief ChronoHydroCoupler: bridges Chrono bodies and HydroForces.
 * 
 * Extracts state from Chrono bodies, evaluates forces via HydroForces,
 * and (in future) applies forces back to Chrono bodies. This is the
 * adapter between Chrono's physics engine and the Chrono-free HydroForces.
 */
class ChronoHydroCoupler {
public:
    /**
     * @brief Constructor.
     * 
     * @param hydro_forces Shared pointer to HydroForces for force evaluation
     * @param bodies Vector of Chrono body pointers
     */
    ChronoHydroCoupler(std::shared_ptr<HydroForces> hydro_forces,
                       std::vector<std::shared_ptr<chrono::ChBody>> bodies);

    /**
     * @brief Evaluate hydrodynamic forces.
     * 
     * Extracts state from Chrono bodies, calls HydroForces to compute
     * forces, and returns the result.
     * 
     * @param time Current simulation time
     * @return BodyForces Total forces (one GeneralizedForce per body, 6 DOF each)
     */
    BodyForces Evaluate(double time);

    /**
     * @brief Apply forces to Chrono bodies.
     * 
     * (Not yet implemented - will be used to push forces back into Chrono)
     * 
     * @param forces Forces to apply
     */
    void ApplyForcesToChrono(const BodyForces& forces);

    /**
     * @brief Get profiling statistics from HydroForces.
     * 
     * Returns cumulative timing and call counts for each component type.
     * 
     * @return HydroForcesProfileStats Profiling statistics
     */
    HydroForcesProfileStats GetProfileStats() const { 
        return hydro_forces_->GetProfileStats(); 
    }

    /**
     * @brief Enable or disable profiling in HydroForces.
     * 
     * @param enabled True to enable profiling, false to disable
     */
    void SetProfilingEnabled(bool enabled) {
        hydro_forces_->SetProfilingEnabled(enabled);
    }

private:
    std::shared_ptr<HydroForces> hydro_forces_;
    std::vector<std::shared_ptr<chrono::ChBody>> bodies_;
};

}  // namespace hydrochrono::hydro

#endif  // HYDROC_COUPLING_CHRONO_COUPLER_H
