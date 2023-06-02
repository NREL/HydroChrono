/*********************************************************************
 * @file  chloadaddedmass.h
 *
 * @brief header file for added mass chload class.
 *********************************************************************/
#pragma once

#include <chrono/core/ChMatrix.h>
#include <chrono/physics/ChBody.h>
#include <hydroc/h5fileinfo.h>
#include <vector>

#include <chrono/physics/ChLoad.h>
#include <chrono/physics/ChSystem.h>

using namespace chrono;
// using namespace chrono::fea;

// =============================================================================
class ChLoadAddedMass : public chrono::ChLoadCustomMultiple {
  public:
    /**
     * @brief Initializes body to have load applied to and added mass matrix from h5 file initialized object.
     *
     * Vector of bodies with hydro forces (and added mass) need to be listed in same order as they are added to the
     * system. Also need to be separate from any bodies without hydro forces (or added mass) and hydro bodies need to be
     * added to system before any bodies without hydro forces applied.
     *
     * @param body_info_struct HydroData::BodyInfo for each body with h5 information including added mass matrix
     * @param bodies vector of Project Chrono bodies to apply added mass to. Must be added to system in same order as in
     * this matrix.
     * @param system pointer to system containing the bodies, used for getting system mass matrix size at any time.
     */
    ChLoadAddedMass(const std::vector<HydroData::BodyInfo>& body_info_struct,
                    std::vector<std::shared_ptr<ChLoadable>>& bodies,
                    ChSystem* system);

    /**
     * @brief "Virtual" copy constructor (covariant return type). Required from chrono inheritance.
     */
    virtual ChLoadAddedMass* Clone() const override { return new ChLoadAddedMass(*this); }

    /**
     * @brief Compute Q, the generalized load.
     *
     * From Chrono documentation:
     * In this case, it computes the quadratic (centrifugal, gyroscopic) terms.
     * Signs are negative as Q assumed at right hand side, so Q= -Fgyro -Fcentrifugal
     * Called automatically at each Update().
     * The M*a term is not added: to this end one could use LoadIntLoadResidual_Mv afterward.
     */
    virtual void ComputeQ(ChState* state_x,      ///< state position to evaluate Q
                          ChStateDelta* state_w  ///< state speed to evaluate Q
                          ) override {}

    /**
     * @brief This is the function that sets the infinite added mass matrix every timestep.
     *
     * From Chrono docs:
     * For efficiency reasons, do not let the parent class do automatic differentiation
     * to compute the R, K matrices. Use analytic expressions instead. For example, R is
     * the well known gyroscopic damping matrix. Also, compute the M matrix.
     *
     * @param state_x state position to evaluate jacobians
     * @param state_w state speed to evaluate jacobians
     * @param mK result -dQ/dx
     * @param mR result -dQ/dv
     * @param mM result -dQ/da
     */
    virtual void ComputeJacobian(ChState* state_x,
                                 ChStateDelta* state_w,
                                 ChMatrixRef mK,
                                 ChMatrixRef mR,
                                 ChMatrixRef mM) override;

    /**
     * @brief Computes LoadIntLoadResidual_Mv for vector w, const c, and vector R. Also carried over from chrono
     * inheritance.
     *
     * Note R here is vector, and is not R gyroscopic damping matrix from ComputeJacobian.
     *  Just for efficiency, override the default LoadIntLoadResidual_Mv, because we can do this in a simplified way.
     *
     * @param R result: the R residual, R += c*M*w
     * @param w the w vector
     * @param c a scaling factor
     */
    virtual void LoadIntLoadResidual_Mv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) override;

  private:
    ChSystem* system;
    ChMatrixDynamic<double> infinite_added_mass;  ///< added mass at infinite frequency in global coordinates
    ChMatrixDynamic<double>
        infinite_added_mass_system;  ///< added mass at infinite frequency in global coordinates (system matrix)
    virtual bool IsStiff() override { return true; }  // this to force the use of the inertial M, R and K matrices
};