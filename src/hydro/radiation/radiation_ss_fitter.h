/*********************************************************************
 * @file  radiation_ss_fitter.h
 * @brief State-space fitter for radiation impulse response kernels.
 *
 * This class implements an algorithm to fit a radiation impulse response
 * K(t) using a state-space realization, matching the WEC-Sim/BEMIO approach.
 *
 * The impulse response is approximated as:
 *     K(t) ≈ C * exp(A*t) * B + D*δ(t)
 *
 * where A, B, C, D are state-space matrices obtained via Hankel-SVD.
 *
 * ALGORITHM OVERVIEW (matches WEC-Sim radiationIRFSS):
 *   1. Build Hankel matrix from scaled kernel samples
 *   2. Perform SVD to get balanced realization
 *   3. Form discrete-time A, B, C, D matrices
 *   4. Convert discrete to continuous using bilinear (Tustin) transform
 *   5. Compute R² to assess fit quality
 *   6. Iterate on order until R² threshold is met or max_order is reached
 *
 * USAGE:
 *   StateSpaceOptions opts;
 *   opts.max_order = 10;
 *   opts.r2_threshold = 0.95;
 *
 *   RadiationStateSpaceFitter fitter(opts);
 *   StateSpaceFitResult result = fitter.FitKernel(K, dt);
 *
 *   // Create model for simulation
 *   auto model = RadiationStateSpaceModel::FromFitResult(result);
 *
 * NOTE: This class is Chrono-free and operates on pure Eigen matrices.
 *********************************************************************/

#ifndef HYDRO_RADIATION_SS_FITTER_H
#define HYDRO_RADIATION_SS_FITTER_H

#include <hydroc/radiation/radiation_types.h>
#include <Eigen/Dense>

namespace hydrochrono::hydro {

/**
 * @brief Result of fitting a state-space model to a kernel.
 *
 * Stores the continuous-time state-space matrices A, B, C, D such that:
 *     K(t) = C * exp(A*t) * B
 *
 * Use RadiationStateSpaceModel::FromFitResult() to convert this to a
 * simulation-ready model with decomposed modes.
 */
struct StateSpaceFitResult {
    Eigen::MatrixXd A;      ///< State matrix [order × order]
    Eigen::VectorXd B;      ///< Input vector [order]
    Eigen::RowVectorXd C;   ///< Output vector [1 × order]
    double D = 0.0;         ///< Feedthrough scalar (usually 0 for radiation)

    double r2 = 0.0;        ///< Coefficient of determination (fit quality, 0-1)
    int order = 0;          ///< State-space order (dimension of A)

    /**
     * @brief Check if the fit is valid (has at least one mode).
     */
    bool IsValid() const { return order > 0; }
};

/**
 * @brief Fitter for state-space approximation of radiation kernels.
 *
 * Given a discrete impulse response kernel K[0], K[1], ..., K[N-1]
 * sampled at intervals dt, this fitter finds a state-space model
 * that approximates the kernel.
 */
class RadiationStateSpaceFitter {
public:
    /**
     * @brief Construct a fitter with given options.
     *
     * @param opts Options controlling max_order and r2_threshold
     */
    explicit RadiationStateSpaceFitter(const StateSpaceOptions& opts);

    /**
     * @brief Fit a single SISO kernel to a state-space model.
     *
     * The kernel K should be the impulse response for a single DOF-pair,
     * sampled at uniform time intervals dt.
     *
     * @param K Kernel samples [N], where K[k] is the value at time k*dt
     * @param dt Time step between samples (seconds)
     *
     * @return StateSpaceFitResult containing A, B, C, D matrices
     *
     * @note Returns order=0 if the kernel is too short or essentially zero.
     */
    StateSpaceFitResult FitKernel(const Eigen::VectorXd& K, double dt) const;

    /**
     * @brief Reconstruct a kernel from state-space fit result.
     *
     * Computes K(t) = C * exp(A*t) * B for each time sample.
     *
     * @param result Fit result containing A, B, C matrices
     * @param dt Time step
     * @param num_samples Number of samples to generate
     *
     * @return Reconstructed kernel [num_samples]
     */
    static Eigen::VectorXd ReconstructKernel(
        const StateSpaceFitResult& result,
        double dt,
        int num_samples);

    /**
     * @brief Compute R² (coefficient of determination) for a fit.
     *
     * R² = 1 - SS_res / SS_tot
     *
     * where SS_res = Σ(K - K_fit)² and SS_tot = Σ(K - mean(K))²
     *
     * @param K_actual Original kernel
     * @param K_fit Reconstructed kernel
     *
     * @return R² value (1.0 = perfect fit, 0.0 = no better than mean)
     */
    static double ComputeR2(const Eigen::VectorXd& K_actual,
                            const Eigen::VectorXd& K_fit);

private:
    StateSpaceOptions opts_;

    /**
     * @brief Internal: Fit with pre-computed SVD for efficiency.
     *
     * This method reuses a single SVD decomposition across multiple order
     * attempts, significantly reducing computation time.
     */
    StateSpaceFitResult FitFromSVD(const Eigen::VectorXd& K,
                                   double dt,
                                   int order,
                                   const Eigen::VectorXd& y,
                                   const Eigen::MatrixXd& U,
                                   const Eigen::MatrixXd& V,
                                   const Eigen::VectorXd& svh,
                                   int hankel_size) const;
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_RADIATION_SS_FITTER_H
