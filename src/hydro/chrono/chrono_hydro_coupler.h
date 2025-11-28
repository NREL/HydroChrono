/*********************************************************************
 * @file  chrono_hydro_coupler.h
 * @brief ChronoHydroCoupler: bridges Chrono bodies and HydroSystem.
 *********************************************************************/

#ifndef HYDRO_CHRONO_CHRONO_HYDRO_COUPLER_H
#define HYDRO_CHRONO_CHRONO_HYDRO_COUPLER_H

#include <hydroc/system_state.h>
#include "../core/hydro_system.h"
#include "../core/chrono_state_utils.h"
#include <chrono/physics/ChBody.h>
#include <vector>
#include <memory>

namespace hydrochrono::hydro {

/**
 * @brief ChronoHydroCoupler: bridges Chrono bodies and HydroSystem.
 * 
 * Extracts state from Chrono bodies, evaluates forces via HydroSystem,
 * and (in future) applies forces back to Chrono bodies. This is the
 * adapter between Chrono's physics engine and the Chrono-free HydroSystem.
 */
class ChronoHydroCoupler {
public:
    /**
     * @brief Constructor.
     * 
     * @param hydro_system Shared pointer to HydroSystem for force evaluation
     * @param bodies Vector of Chrono body pointers
     */
    ChronoHydroCoupler(std::shared_ptr<HydroSystem> hydro_system,
                       std::vector<std::shared_ptr<chrono::ChBody>> bodies);

    /**
     * @brief Evaluate hydrodynamic forces.
     * 
     * Extracts state from Chrono bodies, calls HydroSystem to compute
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
     * @brief Get profiling statistics from HydroSystem.
     * 
     * Returns cumulative timing and call counts for each component type.
     * 
     * @return HydroSystemProfileStats Profiling statistics
     */
    HydroSystemProfileStats GetProfileStats() const { 
        return hydro_system_->GetProfileStats(); 
    }

    /**
     * @brief Enable or disable profiling in HydroSystem.
     * 
     * @param enabled True to enable profiling, false to disable
     */
    void SetProfilingEnabled(bool enabled) {
        hydro_system_->SetProfilingEnabled(enabled);
    }

private:
    std::shared_ptr<HydroSystem> hydro_system_;
    std::vector<std::shared_ptr<chrono::ChBody>> bodies_;
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_CHRONO_CHRONO_HYDRO_COUPLER_H

