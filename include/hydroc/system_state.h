/*********************************************************************
 * @file  system_state.h
 * @brief Common structs for body and system hydrodynamic state.
 *********************************************************************/

#ifndef HYDROC_SYSTEM_STATE_H
#define HYDROC_SYSTEM_STATE_H

#include <Eigen/Dense>
#include <vector>
#include <hydroc/hydro_types.h>

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

#endif  // HYDROC_SYSTEM_STATE_H

