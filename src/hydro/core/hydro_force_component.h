/*********************************************************************
 * @file  hydro_force_component.h
 * @brief Interface for hydrodynamic force components.
 *********************************************************************/

#ifndef HYDRO_CORE_HYDRO_FORCE_COMPONENT_H
#define HYDRO_CORE_HYDRO_FORCE_COMPONENT_H

#include <hydroc/system_state.h>

namespace hydrochrono::hydro {

/**
 * @brief Interface for computing hydrodynamic forces.
 * 
 * All force components (hydrostatics, radiation, excitation) implement this
 * interface to provide a consistent way to compute forces.
 * 
 * @note Components may maintain internal state (e.g., velocity history for
 *       radiation damping), but Compute() should be side-effect-free
 *       outside the component instance.
 */
class IHydroForceComponent {
public:
    virtual ~IHydroForceComponent() = default;

    /**
     * @brief Compute force contribution and add to inout_forces.
     * 
     * Computes the force contribution for the given system state and time,
     * adding the result to the provided forces vector.
     * 
     * @param state Current system state (positions, velocities for all bodies)
     * @param time Current simulation time
     * @param inout_forces Force vector to add contribution to (one GeneralizedForce per body)
     */
    virtual void Compute(const SystemState& state,
                        double time,
                        BodyForces& inout_forces) = 0;
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_CORE_HYDRO_FORCE_COMPONENT_H

