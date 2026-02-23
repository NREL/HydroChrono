/*********************************************************************
 * @file  force_component.h
 * @brief Interface for hydrodynamic force components.
 *
 * MAIN TYPES:
 *   - IHydroForceComponent: Abstract interface for force computation
 *   - HydroComponentType: Enum for profiling (Hydrostatics/Radiation/Excitation)
 *
 * ROLE: All force components implement IHydroForceComponent. HydroSystem
 * owns a collection of components and calls Compute() on each.
 *********************************************************************/

#ifndef HYDROC_CORE_FORCE_COMPONENT_H
#define HYDROC_CORE_FORCE_COMPONENT_H

#include <hydroc/core/system_state.h>

namespace hydrochrono::hydro {

/**
 * @brief Identifies the type of a hydrodynamic force component.
 * 
 * Used for profiling to categorize execution time by component type.
 */
enum class HydroComponentType {
    Hydrostatics,   ///< Hydrostatic restoring forces and buoyancy
    Radiation,      ///< Radiation damping (RIRF convolution)
    Excitation,     ///< Wave excitation forces
    Mooring         ///< Mooring forces (e.g. MoorDyn)
};

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
     * @brief Get the type of this component.
     * 
     * Used for profiling to categorize execution time by component type.
     * 
     * @return Component type (Hydrostatics, Radiation, or Excitation)
     */
    virtual HydroComponentType Type() const = 0;

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

#endif  // HYDROC_CORE_FORCE_COMPONENT_H

