/*********************************************************************
 * @file  hydro_system.h
 * @brief HydroSystem façade: orchestrates force components to compute total forces.
 *********************************************************************/

#ifndef HYDRO_CORE_HYDRO_SYSTEM_H
#define HYDRO_CORE_HYDRO_SYSTEM_H

#include "system_state.h"
#include "hydro_force_component.h"
#include <vector>
#include <memory>

namespace hydrochrono::hydro {

/**
 * @brief HydroSystem: orchestrates force components to compute total forces.
 * 
 * Owns a collection of force components (hydrostatics, radiation, excitation)
 * and evaluates their combined contributions given system state and time.
 * This is a Chrono-free façade for force computation.
 */
class HydroSystem {
public:
    /**
     * @brief Constructor.
     * 
     * @param num_bodies Number of bodies in the system
     * @param components Vector of force components (ownership transferred via move)
     */
    HydroSystem(int num_bodies,
                std::vector<std::unique_ptr<IHydroForceComponent>> components);

    /**
     * @brief Evaluate total hydrodynamic forces.
     * 
     * Computes force contributions from all components for the given
     * system state and time, accumulating them into a single result.
     * 
     * @param state Current system state (positions, velocities for all bodies)
     * @param time Current simulation time
     * @return BodyForces Total forces (one GeneralizedForce per body, 6 DOF each)
     */
    BodyForces Evaluate(const SystemState& state, double time) const;

    /**
     * @brief Get number of bodies in the system.
     * 
     * @return Number of bodies
     */
    int num_bodies() const { return num_bodies_; }

private:
    int num_bodies_;
    std::vector<std::unique_ptr<IHydroForceComponent>> components_;
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_CORE_HYDRO_SYSTEM_H

