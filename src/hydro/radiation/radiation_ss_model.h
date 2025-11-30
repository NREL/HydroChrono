/*********************************************************************
 * @file  radiation_ss_model.h
 * @brief State-space model for radiation damping with full mode support.
 *
 * This class implements a state-space approximation of the radiation
 * impulse response function (RIRF) supporting both:
 *   - Pure exponential modes: K(t) = H * exp(-α * t)
 *   - Oscillatory modes: K(t) = exp(-α*t) * (H_c*cos(ω*t) + H_s*sin(ω*t))
 *
 * This enables O(1) per-timestep computation for any kernel, including those
 * with damped oscillations (complex conjugate poles), which are typical for
 * hydrodynamic radiation problems.
 *
 * THEORY:
 *   Pure exponential mode (real eigenvalue λ = -α):
 *       ż(t) = -α z(t) + v(t)
 *       f_contribution = H * z
 *
 *   Oscillatory mode (complex conjugate pair λ = -α ± jω):
 *       d/dt [z_c]   = [-α  -ω] [z_c]   + [b_c] * v(t)
 *       d/dt [z_s]     [ ω  -α] [z_s]     [b_s]
 *       f_contribution = H_c * z_c + H_s * z_s
 *
 * Both use exact exponential integration for unconditional stability.
 *
 * USAGE:
 *   // From fit result (automatic mode extraction)
 *   auto model = RadiationStateSpaceModel::FromFitResult(result);
 *   model.Reset();
 *   model.Step(dt, v);
 *   double f = model.GetForce();
 *
 *   // Or manually construct modes
 *   RadiationStateSpaceModel model;
 *   model.AddOscillatoryMode(0.8, 3.0, 1.0, 0.0, 5.0, 0.0);
 *   model.Reset();
 *
 * NOTE: This class is Chrono-free and operates on pure Eigen/scalars.
 *       For multi-DOF systems, use one model per (i,j) DOF pair.
 *********************************************************************/

#ifndef HYDRO_RADIATION_SS_MODEL_H
#define HYDRO_RADIATION_SS_MODEL_H

#include <Eigen/Dense>
#include <vector>

namespace hydrochrono::hydro {

// Forward declaration
struct StateSpaceFitResult;

/**
 * @brief A single pure exponential mode (real eigenvalue).
 *
 * ODE: z' = -α*z + b*v
 * Force contribution: H * z
 * Impulse response: K(t) = H * b * exp(-α*t)
 */
struct ExponentialMode {
    double alpha;   ///< Decay rate > 0
    double b;       ///< Input gain (how velocity affects state)
    double H;       ///< Output gain (how state affects force)
    double z;       ///< Internal state
    
    ExponentialMode() : alpha(1.0), b(1.0), H(1.0), z(0.0) {}
    ExponentialMode(double a, double b_, double h) : alpha(a), b(b_), H(h), z(0.0) {}
};

/**
 * @brief An oscillatory mode (complex conjugate eigenvalue pair).
 *
 * ODE: d/dt [z_c] = [-α  -ω] [z_c] + [b_c] * v
 *      d/dt [z_s]   [ ω  -α] [z_s]   [b_s]
 * Force contribution: H_c * z_c + H_s * z_s
 * Impulse response: K(t) = exp(-α*t) * (H_c*cos(ω*t) + H_s*sin(ω*t))
 */
struct OscillatoryMode {
    double alpha;   ///< Decay rate > 0
    double omega;   ///< Oscillation frequency > 0
    double b_c;     ///< Input gain for cosine component
    double b_s;     ///< Input gain for sine component
    double H_c;     ///< Output gain for cosine component
    double H_s;     ///< Output gain for sine component
    double z_c;     ///< Internal state (cosine component)
    double z_s;     ///< Internal state (sine component)
    
    OscillatoryMode() : alpha(1.0), omega(1.0), b_c(1.0), b_s(0.0), 
                        H_c(1.0), H_s(0.0), z_c(0.0), z_s(0.0) {}
};

/**
 * @brief State-space model for radiation damping (SISO).
 *
 * Supports both pure exponential and oscillatory modes for O(1) per-timestep.
 * Use FromFitResult() to construct from a StateSpaceFitResult.
 *
 * For multi-DOF systems (e.g., 6-DOF body), create one RadiationStateSpaceModel
 * per (input_dof, output_dof) pair, then sum the forces.
 */
class RadiationStateSpaceModel {
public:
    /**
     * @brief Default constructor (empty model).
     */
    RadiationStateSpaceModel() = default;

    /**
     * @brief Create model from a state-space fit result.
     *
     * Automatically extracts real and complex eigenvalues from the A matrix
     * and creates the appropriate mode types.
     *
     * @param result Fit result containing A, B, C matrices
     * @return Model with decomposed modes
     */
    static RadiationStateSpaceModel FromFitResult(const StateSpaceFitResult& result);

    /**
     * @brief Add a pure exponential mode.
     *
     * @param alpha Decay rate > 0
     * @param b Input gain
     * @param H Output gain
     */
    void AddExponentialMode(double alpha, double b, double H);

    /**
     * @brief Add an oscillatory mode.
     *
     * @param alpha Decay rate > 0
     * @param omega Oscillation frequency > 0
     * @param b_c Input gain for cos component
     * @param b_s Input gain for sin component
     * @param H_c Output gain for cos component
     * @param H_s Output gain for sin component
     */
    void AddOscillatoryMode(double alpha, double omega, 
                            double b_c, double b_s,
                            double H_c, double H_s);

    /**
     * @brief Reset all internal states to zero.
     */
    void Reset();

    /**
     * @brief Advance the model by one time step.
     *
     * Uses exact exponential integration:
     * - Exponential modes: z(t+dt) = exp(-α*dt)*z(t) + (1-exp(-α*dt))/α * v
     * - Oscillatory modes: uses rotation matrix with decay
     *
     * Both are unconditionally stable for any dt > 0.
     *
     * @param dt Time step (seconds, > 0)
     * @param v Velocity input (scalar)
     */
    void Step(double dt, double v);

    /**
     * @brief Get the current radiation force.
     *
     * @return Sum of all mode contributions: Σ(H*z) for exp, Σ(H_c*z_c + H_s*z_s) for osc
     */
    double GetForce() const;

    /**
     * @brief Reconstruct the impulse response kernel.
     *
     * Useful for verification: K_reconstructed should match K_original if fit is good.
     *
     * @param dt Time step
     * @param num_samples Number of samples to generate
     * @return Reconstructed kernel K[0], K[1], ..., K[num_samples-1]
     */
    Eigen::VectorXd ReconstructKernel(double dt, int num_samples) const;

    // Accessors for diagnostics
    int num_exp_modes() const { return static_cast<int>(exp_modes_.size()); }
    int num_osc_modes() const { return static_cast<int>(osc_modes_.size()); }
    int total_modes() const { return num_exp_modes() + num_osc_modes(); }
    int total_states() const { return num_exp_modes() + 2 * num_osc_modes(); }

    const std::vector<ExponentialMode>& exp_modes() const { return exp_modes_; }
    const std::vector<OscillatoryMode>& osc_modes() const { return osc_modes_; }

private:
    std::vector<ExponentialMode> exp_modes_;   ///< Pure exponential modes
    std::vector<OscillatoryMode> osc_modes_;   ///< Oscillatory modes
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_RADIATION_SS_MODEL_H
