/*********************************************************************
 * @file  hydrostatics_component.h
 * @brief Hydrostatics force component (restoring forces and buoyancy).
 *********************************************************************/

#ifndef HYDRO_FORCE_COMPONENTS_HYDROSTATICS_COMPONENT_H
#define HYDRO_FORCE_COMPONENTS_HYDROSTATICS_COMPONENT_H

#include "../core/hydro_force_component.h"
#include <hydroc/system_state.h>
#include <Eigen/Dense>
#include <vector>
#include <memory>

class HydroData;

namespace hydrochrono::hydro {

/**
 * @brief Hydrostatics force component (restoring stiffness + buoyancy).
 * 
 * Computes hydrostatic restoring forces from linear stiffness matrix
 * and buoyancy forces at equilibrium position.
 */
class HydrostaticsComponent : public IHydroForceComponent {
public:
    /**
     * @brief Constructor.
     * 
     * @param file_info Reference to HydroData containing stiffness matrices, equilibrium data
     * @param num_bodies Number of bodies in the system
     * @param equilibrium Equilibrium positions [6N] (xyz, rpy per body)
     * @param cb_minus_cg Center of buoyancy minus center of gravity [3N] (xyz per body)
     * @param gravitational_acceleration Gravitational acceleration vector (m/s^2)
     */
    HydrostaticsComponent(const HydroData& file_info,
                     int num_bodies,
                     const std::vector<double>& equilibrium,
                     const std::vector<double>& cb_minus_cg,
                     const Eigen::Vector3d& gravitational_acceleration);

    /**
     * @brief Get the component type.
     * @return HydroComponentType::Hydrostatics
     */
    HydroComponentType Type() const override { return HydroComponentType::Hydrostatics; }

    /**
     * @brief Compute hydrostatic force contribution.
     * 
     * Adds hydrostatic restoring forces and buoyancy forces to inout_forces.
     * Forces computed per body from displacement from equilibrium.
     * 
     * @param state Current system state (positions and orientations)
     * @param time Current simulation time (not used for hydrostatics)
     * @param inout_forces Force vector to add contribution to (one GeneralizedForce per body)
     */
    void Compute(const SystemState& state,
                double time,
                BodyForces& inout_forces) override;

private:
    const HydroData& file_info_;
    int num_bodies_;
    std::vector<double> equilibrium_;
    std::vector<double> cb_minus_cg_;
    Eigen::Vector3d gravitational_acceleration_;
    
    static constexpr int kDofPerBody = 6;
    static constexpr int kDofLinOrRot = 3;
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_FORCE_COMPONENTS_HYDROSTATICS_COMPONENT_H

