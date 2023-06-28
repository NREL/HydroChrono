/*********************************************************************
 * @file  chloadaddedmass.h
 *
 * @brief header file for added mass chload class.
 *********************************************************************/
#include <hydroc/chloadaddedmass.h>

#include <utility>

#include "chrono/physics/ChBody.h"

ChLoadAddedMass::ChLoadAddedMass(const std::vector<HydroData::BodyInfo>& user_h5_body_data,
                                 std::vector<std::shared_ptr<ChLoadable>>& bodies,
                                 ChSystem* system)
    : ChLoadCustomMultiple(bodies), system(system) {
    auto nBodies = bodies.size();

    infinite_added_mass.setZero(6 * nBodies, 6 * nBodies);
    for (int i = 0; i < nBodies; i++) {
        infinite_added_mass.block(i * 6, 0, 6, nBodies * 6) = user_h5_body_data[i].inf_added_mass;
    }

    // initialize added mass matrix for whole system
    infinite_added_mass_system = infinite_added_mass;
}

void ChLoadAddedMass::ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
                                      ChStateDelta* state_w,  ///< state speed to evaluate jacobians
                                      ChMatrixRef mK,         ///< result dQ/dx
                                      ChMatrixRef mR,         ///< result dQ/dv
                                      ChMatrixRef mM          ///< result dQ/da
) {
    // The following ensures that the added mass matrix matches the size of the ChSystem mass matrix. It is necessary
    // for systems that have both hydro and non-hydro bodies when adding a system-wide load.
    // @todo if possible, remove hack by using initialiazer function called AFTER the ChSystem is assembled.

    // get mass matrix size
    auto mmrows = system->GetNcoords_w();  // number of rows in system mass matrix
    // check if ChSystem mass matrix size different from added mass matrix size
    if (mmrows != infinite_added_mass_system.rows() && mmrows > 0) {
        // initialize/update system matrix;
        infinite_added_mass_system.setZero(mmrows, mmrows);
        auto amrows                                            = infinite_added_mass.rows();
        infinite_added_mass_system.block(0, 0, amrows, amrows) = infinite_added_mass;
    }
    // set mass matrix here
    jacobians->M = infinite_added_mass_system;

    // R gyroscopic damping matrix terms (6Nx6N)
    // 0 for added mass
    jacobians->R.setZero();

    // K inertial stiffness matrix terms (6Nx6N)
    // 0 for added mass
    jacobians->K.setZero();
}

void ChLoadAddedMass::LoadIntLoadResidual_Mv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) {
    if (!this->jacobians) return;

    // if (!loadable->IsSubBlockActive(0))
    //	return;

    // R+=c*M*a
    // segment gives the chunk of vector starting at the first argument, and going for as many elements as the second
    // argument... in this case, segment gets the 3vector starting at the 0th DOF's offset (ie 0)
    // R.segment(loadable->GetSubBlockOffset(0), 3) += c * (this->mass * (a_x + chrono::Vcross(a_w,
    // this->c_m))).eigen();
    // in this case, segment gets the 3vector starting at the 0th DOF's + 3 offset (ie 3)
    // R.segment(loadable->GetSubBlockOffset(0) + 3, 3) += c * (this->mass * chrono::Vcross(this->c_m, a_x) + this->I *
    // a_w).eigen();
    // since R is a vector, we can probably just do R += C*M*a with no need to separate w into a_x and a_w above
    R += c * jacobians->M * w;
}