#include <hydroc/chloadaddedmass.h>

#include "chrono/physics/ChBody.h"
#include <utility>


// =============================================================================
// ChLoadAddedMass Class Definitions
// =============================================================================



/*******************************************************************************
* ChLoadAddedMass constructor
* initializes body to have load applied to and added mass matrix from h5 file object
*******************************************************************************/
ChLoadAddedMass::ChLoadAddedMass(const std::vector<H5FileInfo>& user_h5_body_data,
								 std::vector<std::shared_ptr<ChLoadable>>& bodies)
								     : ChLoadCustomMultiple(bodies) { 
	nBodies = bodies.size();

    infinite_added_mass.setZero(6 * nBodies, 6 * nBodies);
	for (int i = 0; i < nBodies; i++) {
		infinite_added_mass.block(i * 6, 0, 6, nBodies * 6) = user_h5_body_data[i].GetInfAddedMassMatrix();
	}
}

/*******************************************************************************
* ChLoadAddedMass::ComputeJacobian()
* Computes Jacobian for load, in this case just the mass matrix is initialized
* as the added mass matrix
*******************************************************************************/
void ChLoadAddedMass::ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
	ChStateDelta* state_w,  ///< state speed to evaluate jacobians
	ChMatrixRef mK,         ///< result dQ/dx
	ChMatrixRef mR,         ///< result dQ/dv
	ChMatrixRef mM          ///< result dQ/da
) {
	//set mass matrix here
	jacobians->M = infinite_added_mass;

	// R gyroscopic damping matrix terms (6Nx6N)
	// 0 for added mass
	jacobians->R.setZero();

	// K inertial stiffness matrix terms (6Nx6N)
	// 0 for added mass
	jacobians->K.setZero();
}

/*******************************************************************************
* ChLoadAddedMass::LoadIntLoadResidual_Mv()
* Computes LoadIntLoadResidual_Mv for vector w, const c, and vector R
* Note R here is vector, and is not R gyroscopic damping matrix from ComputeJacobian
*******************************************************************************/
void ChLoadAddedMass::LoadIntLoadResidual_Mv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) {
	if (!this->jacobians)
		return;

	//if (!loadable->IsSubBlockActive(0))
	//	return;

	// R+=c*M*a
	// segment gives the chunk of vector starting at the first argument, and going for as many elements as the second argument...
	// in this case, segment gets the 3vector starting at the 0th DOF's offset (ie 0)
	//R.segment(loadable->GetSubBlockOffset(0), 3) += c * (this->mass * (a_x + chrono::Vcross(a_w, this->c_m))).eigen();
	// in this case, segment gets the 3vector starting at the 0th DOF's + 3 offset (ie 3)
	//R.segment(loadable->GetSubBlockOffset(0) + 3, 3) += c * (this->mass * chrono::Vcross(this->c_m, a_x) + this->I * a_w).eigen();
	// since R is a vector, we can probably just do R += C*M*a with no need to separate w into a_x and a_w above
	R += c * jacobians->M * w;
}