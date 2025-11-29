/*********************************************************************
 * @file  system_state.h
 * @brief Chrono-free state structs for hydrodynamic computation.
 *
 * MAIN TYPES:
 *   - BodyState: Position, orientation, velocities for one body
 *   - SystemState: Vector of BodyState for all bodies
 *
 * ROLE: Input to HydroSystem and force components. Decouples core
 * hydrodynamics from Chrono's ChBody representation.
 *********************************************************************/

#ifndef HYDROC_CORE_SYSTEM_STATE_H
#define HYDROC_CORE_SYSTEM_STATE_H

#include <Eigen/Dense>
#include <vector>
#include <hydroc/core/hydro_types.h>

namespace hydrochrono::hydro {

/**
 * @brief State of a single body (pose and velocities).
 * 
 * Contains position, orientation (as RPY angles), linear velocity,
 * and angular velocity for one body in the system.
 */
struct BodyState {
    Eigen::Vector3d position;        // Position in world frame (x, y, z)
    Eigen::Vector3d orientation_rpy; // Orientation as roll-pitch-yaw angles (radians)
    Eigen::Vector3d linear_velocity;  // Linear velocity in world frame (m/s)
    Eigen::Vector3d angular_velocity; // Angular velocity in parent frame (rad/s)
};

/**
 * @brief State of all bodies in the hydrodynamic system.
 * 
 * Contains a vector of BodyState, one per body.
 */
struct SystemState {
    std::vector<BodyState> bodies;
};

}  // namespace hydrochrono::hydro

#endif  // HYDROC_CORE_SYSTEM_STATE_H

