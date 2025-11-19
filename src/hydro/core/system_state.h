/*********************************************************************
 * @file  system_state.h
 * @brief Common structs for body and system hydrodynamic state.
 *********************************************************************/

#ifndef HYDRO_CORE_SYSTEM_STATE_H
#define HYDRO_CORE_SYSTEM_STATE_H

#include <Eigen/Dense>
#include <vector>

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

/**
 * @brief Generalized force vector for one body (6 DOF).
 * 
 * Force vector [Fx, Fy, Fz, Mx, My, Mz] for surge, sway, heave,
 * roll, pitch, yaw degrees of freedom.
 */
using GeneralizedForce = Eigen::VectorXd;

/**
 * @brief Forces for all bodies in the system.
 * 
 * Vector of GeneralizedForce, one per body.
 */
using BodyForces = std::vector<GeneralizedForce>;

}  // namespace hydrochrono::hydro

#endif  // HYDRO_CORE_SYSTEM_STATE_H

