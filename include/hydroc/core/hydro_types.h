/*********************************************************************
 * @file  hydro_types.h
 * @brief Core type aliases for hydrodynamic forces.
 *
 * MAIN TYPES:
 *   - GeneralizedForce: 6-DOF force vector for one body (Eigen::VectorXd)
 *   - BodyForces: Vector of GeneralizedForce for all bodies
 *
 * ROLE: Canonical definitions used throughout the hydrodynamics core.
 *********************************************************************/

#ifndef HYDROC_CORE_HYDRO_TYPES_H
#define HYDROC_CORE_HYDRO_TYPES_H

#include <Eigen/Dense>
#include <vector>

namespace hydrochrono::hydro {

/**
 * @brief Generalized force vector for one body (6 DOF).
 * 
 * Force vector [Fx, Fy, Fz, Mx, My, Mz] for surge, sway, heave,
 * roll, pitch, yaw degrees of freedom.
 */
using GeneralizedForce = Eigen::VectorXd;   //// RADU - why is this not an Eigen vector of fixed size 6, or better yet 2 vectors?

/**
 * @brief Forces for all bodies in the system.
 * 
 * Vector of GeneralizedForce, one per body.
 */
using BodyForces = std::vector<GeneralizedForce>;

}  // namespace hydrochrono::hydro

#endif  // HYDROC_CORE_HYDRO_TYPES_H

