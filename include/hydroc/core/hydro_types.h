/*********************************************************************
 * @file  hydro_types.h
 * @brief Canonical type definitions for hydrodynamic force types.
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
using GeneralizedForce = Eigen::VectorXd;

/**
 * @brief Forces for all bodies in the system.
 * 
 * Vector of GeneralizedForce, one per body.
 */
using BodyForces = std::vector<GeneralizedForce>;

}  // namespace hydrochrono::hydro

#endif  // HYDROC_CORE_HYDRO_TYPES_H

