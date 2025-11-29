/*********************************************************************
 * @file  hydro_system.h
 * @brief HydroSystem façade: orchestrates force components to compute total forces.
 *********************************************************************/

#ifndef HYDROC_CORE_HYDRO_SYSTEM_H
#define HYDROC_CORE_HYDRO_SYSTEM_H

#include <hydroc/core/system_state.h>
#include <hydroc/core/force_component.h>
#include <vector>
#include <memory>

namespace hydrochrono::hydro {

/**
 * @brief Profiling statistics for HydroSystem components.
 * 
 * Tracks cumulative execution time and call counts for each component type.
 */
struct HydroSystemProfileStats {
    double hydrostatics_seconds = 0.0;
    double radiation_seconds    = 0.0;
    double excitation_seconds   = 0.0;
    int hydrostatics_calls      = 0;
    int radiation_calls         = 0;
    int excitation_calls        = 0;
    
    /// Reset all counters to zero
    void Reset() {
        hydrostatics_seconds = radiation_seconds = excitation_seconds = 0.0;
        hydrostatics_calls = radiation_calls = excitation_calls = 0;
    }
};

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
     * Also updates profiling statistics for each component type.
     * 
     * @param state Current system state (positions, velocities for all bodies)
     * @param time Current simulation time
     * @return BodyForces Total forces (one GeneralizedForce per body, 6 DOF each)
     */
    BodyForces Evaluate(const SystemState& state, double time);

    /**
     * @brief Get number of bodies in the system.
     * 
     * @return Number of bodies
     */
    int num_bodies() const { return num_bodies_; }

    /**
     * @brief Get profiling statistics.
     * 
     * Returns cumulative timing and call counts for each component type.
     * 
     * @return HydroSystemProfileStats Profiling statistics
     */
    HydroSystemProfileStats GetProfileStats() const { return profile_stats_; }

    /**
     * @brief Reset profiling statistics.
     * 
     * Clears all timing and call counts.
     */
    void ResetProfileStats() { profile_stats_.Reset(); }

    /**
     * @brief Enable or disable profiling.
     * 
     * When disabled (default), Evaluate() skips timing measurements.
     * 
     * @param enabled True to enable profiling, false to disable
     */
    void SetProfilingEnabled(bool enabled) { profiling_enabled_ = enabled; }

    /**
     * @brief Check if profiling is enabled.
     * 
     * @return True if profiling is enabled
     */
    bool ProfilingEnabled() const { return profiling_enabled_; }

private:
    int num_bodies_;
    std::vector<std::unique_ptr<IHydroForceComponent>> components_;
    mutable HydroSystemProfileStats profile_stats_;
    bool profiling_enabled_ = false;  ///< Profiling disabled by default
};

}  // namespace hydrochrono::hydro

#endif  // HYDROC_CORE_HYDRO_SYSTEM_H

