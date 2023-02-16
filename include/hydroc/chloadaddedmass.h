#pragma once

#include <vector>
#include <hydroc/h5fileinfo.h>
#include <chrono/core/ChMatrix.h>
#include <chrono/physics/ChBody.h>

#include <chrono/physics/ChLoad.h>

using namespace chrono;
// using namespace chrono::fea;

// =============================================================================
class ChLoadAddedMass : public chrono::ChLoadCustomMultiple {
  public:
    ChLoadAddedMass(const std::vector<H5FileInfo>& file,              ///< h5 file to initialize added mass with
                    std::vector<std::shared_ptr<ChLoadable>>& bodies  ///< objects to apply additional inertia to
    );
    //	/// "Virtual" copy constructor (covariant return type).
    virtual ChLoadAddedMass* Clone() const override { return new ChLoadAddedMass(*this); }
    //
    //	/// Compute Q, the generalized load.
    //	/// In this case, it computes the quadratic (centrifugal, gyroscopic) terms.
    //	/// Signs are negative as Q assumed at right hand side, so Q= -Fgyro -Fcentrifugal
    //	/// Called automatically at each Update().
    //	/// The M*a term is not added: to this end one could use LoadIntLoadResidual_Mv afterward.
    virtual void ComputeQ(ChState* state_x,      ///< state position to evaluate Q
                          ChStateDelta* state_w  ///< state speed to evaluate Q
                          ) override {}

    /// For efficiency reasons, do not let the parent class do automatic differentiation
    /// to compute the R, K matrices. Use analytic expressions instead. For example, R is
    /// the well known gyroscopic damping matrix. Also, compute the M matrix.
    virtual void ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
                                 ChStateDelta* state_w,  ///< state speed to evaluate jacobians
                                 ChMatrixRef mK,         ///< result -dQ/dx
                                 ChMatrixRef mR,         ///< result -dQ/dv
                                 ChMatrixRef mM          ///< result -dQ/da
                                 ) override;

    /// Just for efficiency, override the default LoadIntLoadResidual_Mv, because we can do this in a simplified way.
    virtual void LoadIntLoadResidual_Mv(ChVectorDynamic<>& R,        ///< result: the R residual, R += c*M*w
                                        const ChVectorDynamic<>& w,  ///< the w vector
                                        const double c) override;    ///< a scaling factor

  private:
    ChMatrixDynamic<double> infinite_added_mass;  ///< added mass at infinite frequency in global coordinates
    int nBodies;
    virtual bool IsStiff() override { return true; }  // this to force the use of the inertial M, R and K matrices
};