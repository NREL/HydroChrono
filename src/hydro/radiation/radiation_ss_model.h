/*********************************************************************
 * @file  radiation_ss_model.h
 * @brief State-space model for radiation damping (exponential kernel approximation).
 *
 * This class implements the mathematical core of state-space radiation:
 *   - Approximates the radiation impulse response as a sum of exponentials
 *   - Integrates internal states using an exact exponential scheme
 *   - Returns radiation forces given velocity inputs
 *
 * THEORY:
 *   The radiation kernel K(τ) is approximated as:
 *       K(τ) ≈ Σₘ Hₘ exp(−αₘ τ)
 *
 *   This leads to ODEs for internal states zₘ:
 *       żₘ(t) = v(t) − αₘ zₘ(t)
 *
 *   And radiation forces:
 *       f_rad(t) = Σₘ Hₘ ⊙ zₘ(t)
 *
 *   where ⊙ denotes element-wise multiplication.
 *
 * USAGE:
 *   1. Construct with number of DOFs and modes
 *   2. Set parameters for each mode via SetModeParameters()
 *   3. Call Reset() to zero states
 *   4. In simulation loop: Step(dt, v), then GetForces()
 *
 * NOTE: This class is Chrono-free and operates on pure Eigen vectors.
 *********************************************************************/

#ifndef HYDRO_RADIATION_SS_MODEL_H
#define HYDRO_RADIATION_SS_MODEL_H

#include <Eigen/Dense>

namespace hydrochrono::hydro {

/**
 * @brief State-space model for radiation damping.
 *
 * Models radiation forces using a sum of exponential decay modes.
 * Each mode m has:
 *   - A decay rate αₘ > 0
 *   - A gain vector Hₘ of size num_dofs
 *   - An internal state vector zₘ of size num_dofs
 *
 * The model uses an exact exponential integrator for unconditional stability.
 */
class RadiationStateSpaceModel {
public:
    /**
     * @brief Construct a state-space model.
     *
     * @param num_dofs Number of degrees of freedom (e.g., 6 for single body, 12 for two bodies)
     * @param num_modes Number of exponential modes in the approximation
     *
     * @throws std::invalid_argument if num_dofs <= 0 or num_modes <= 0
     */
    RadiationStateSpaceModel(int num_dofs, int num_modes);

    /**
     * @brief Set parameters for a single mode.
     *
     * Must be called for each mode (0 to num_modes-1) before stepping.
     *
     * @param mode_index Index of the mode (0-based)
     * @param alpha Decay rate (must be > 0 for stability)
     * @param h_column Gain vector mapping this mode's state to forces (size num_dofs)
     *
     * @throws std::out_of_range if mode_index is invalid
     * @throws std::invalid_argument if alpha <= 0 or h_column size mismatches
     */
    void SetModeParameters(int mode_index, double alpha, const Eigen::VectorXd& h_column);

    /**
     * @brief Reset all internal states to zero.
     *
     * Call this at the start of a simulation or to reinitialize.
     */
    void Reset();

    /**
     * @brief Advance the model by one time step.
     *
     * Uses exact exponential integration:
     *   zₘ(t+dt) = exp(−αₘ dt) zₘ(t) + [1 − exp(−αₘ dt)] / αₘ · v
     *
     * This scheme is unconditionally stable for any dt > 0 and αₘ > 0.
     *
     * @param dt Time step (seconds, must be > 0)
     * @param v Velocity vector for all DOFs (size num_dofs)
     *
     * @throws std::invalid_argument if dt <= 0 or v size mismatches
     */
    void Step(double dt, const Eigen::VectorXd& v);

    /**
     * @brief Get the current radiation force vector.
     *
     * Computes: f_rad = Σₘ Hₘ ⊙ zₘ
     * where ⊙ is element-wise multiplication.
     *
     * @return Radiation force vector (size num_dofs)
     */
    Eigen::VectorXd GetForces() const;

    // Accessors for testing and diagnostics
    int num_dofs() const { return num_dofs_; }
    int num_modes() const { return num_modes_; }
    const Eigen::VectorXd& alphas() const { return alphas_; }
    const Eigen::MatrixXd& H() const { return H_; }
    const Eigen::MatrixXd& z() const { return z_; }

private:
    int num_dofs_;          ///< Number of degrees of freedom
    int num_modes_;         ///< Number of exponential modes

    Eigen::VectorXd alphas_;  ///< Decay rates [num_modes]
    Eigen::MatrixXd H_;       ///< Gain matrix [num_dofs × num_modes]
    Eigen::MatrixXd z_;       ///< Internal states [num_dofs × num_modes]
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_RADIATION_SS_MODEL_H

