/*********************************************************************
 * @file  added_mass.h
 * @brief ChLoadAddedMass: Chrono load for infinite-frequency added mass.
 *
 * STATUS: LEGACY IMPLEMENTATION -- retained for reference and as a
 * switchable alternative. The default added-mass path now uses
 * Chrono's built-in ChLoadHydrodynamics.
 *
 * To activate this implementation, define HYDROCHRONO_USE_LEGACY_ADDED_MASS
 * in hydro_system.h (see toggle comment near the top of that file).
 *
 * MAIN TYPES:
 *   - ChLoadAddedMass: ChLoadCustomMultiple subclass that applies
 *     added mass from HDF5 hydrodynamic data to Chrono bodies.
 *
 * ROLE: Chrono coupling layer. Provides a ChLoadCustomMultiple-based
 * approach to applying added mass forces to the Chrono simulation system.
 *********************************************************************/

#ifndef HYDROC_COUPLING_ADDED_MASS_H
#define HYDROC_COUPLING_ADDED_MASS_H

#include <chrono/core/ChMatrix.h>
#include <chrono/physics/ChBody.h>
#include <hydroc/io/h5_reader.h>
#include <vector>

#include <chrono/physics/ChLoad.h>
#include <chrono/physics/ChSystem.h>

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
                    std::vector<std::shared_ptr<chrono::ChLoadable>>& bodies,
                    chrono::ChSystem* system);

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
    virtual void ComputeQ(chrono::ChState* state_x,      ///< state position to evaluate Q
                          chrono::ChStateDelta* state_w  ///< state speed to evaluate Q
                          ) override {}

    /**
     * @brief This is the function that sets the infinite added mass matrix every timestep.
     *
     * From Chrono docs:
     * Compute the K=-dQ/dx, R=-dQ/dv, M=-dQ/da Jacobians.
     * Implementation in a derived class should load the Jacobian matrices K, R, M in the structure 'm_jacobians'.
     * Note the sign that is flipped because we assume equations are written with Q moved to left-hand side.
     *
     * @param state_x state position to evaluate jacobians
     * @param state_w state speed to evaluate jacobians
     */
    virtual void ComputeJacobian(chrono::ChState* state_x, chrono::ChStateDelta* state_w) override;

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
    virtual void LoadIntLoadResidual_Mv(chrono::ChVectorDynamic<>& R,
                                        const chrono::ChVectorDynamic<>& w,
                                        const double c) override;

  private:
    chrono::ChSystem* system;
    chrono::ChMatrixDynamic<double> infinite_added_mass;  ///< added mass at infinite frequency in global coordinates
    chrono::ChMatrixDynamic<double>
        infinite_added_mass_system;  ///< added mass at infinite frequency in global coordinates (system matrix)
    virtual bool IsStiff() override { return true; }  // this to force the use of the inertial M, R and K matrices
};

#endif  // HYDROC_COUPLING_ADDED_MASS_H

