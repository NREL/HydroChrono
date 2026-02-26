/*********************************************************************
 * @file  radiation_types.h
 * @brief Public types for radiation damping configuration.
 *
 * MAIN TYPES:
 *   - RadiationMethod: Selects the radiation damping calculation method
 *   - RadiationKernelProcessing: Composable smoothing + tapering for RIRF kernels
 *   - StateSpaceOptions: Settings for state-space radiation approximation
 *********************************************************************/

#ifndef HYDROC_RADIATION_TYPES_H
#define HYDROC_RADIATION_TYPES_H

#include <string>

namespace hydrochrono::hydro {

/**
 * @brief Radiation damping calculation method.
 * 
 * Selects the overall approach for computing radiation damping forces:
 *   - kRirfConvolution: Direct convolution of RIRF kernels with velocity history (default)
 *   - kStateSpace: State-space approximation using exponential decay modes
 * 
 * RIRF truncation (SetRadiationTruncationTime) is applied before either method.
 * Kernel processing (smoothing/tapering) is applied only for kRirfConvolution.
 */
enum class RadiationMethod {
    kRirfConvolution,  ///< Direct RIRF convolution (default, existing method)
    kStateSpace        ///< State-space approximation (exponential decay modes)
};

/**
 * @brief Composable RIRF kernel processing options (smoothing and/or tapering).
 * 
 * Applied after RIRF truncation (which is handled separately via
 * HydroSystem::SetRadiationTruncationTime). Smoothing is applied first,
 * then tapering.
 * 
 * Defaults are no-ops: smoothing is "none" and taper is disabled. If neither
 * is configured, no kernel processing occurs.
 */
struct RadiationKernelProcessing {
    /// Smoothing type: "none" (default), "sg" (Savitzky-Golay), or "moving_average"
    std::string smoothing_type = "none";
    int smoothing_window = 5;              ///< Smoothing window length (odd, >= 3)

    /// Tapering: when enabled, applies a half-cosine taper window near the end of the RIRF
    bool taper_enabled = false;
    double taper_start_percent = 0.8;      ///< Start taper at this fraction of RIRF length
    double taper_end_percent = 1.0;        ///< End taper at this fraction of RIRF length
    double taper_final_amplitude = 0.0;    ///< Final amplitude (0.0 = zero, 1.0 = no change)

    bool export_csv = false;               ///< Export before/after CSV for diagnostics

    /// Returns true if any processing (smoothing or tapering) is configured.
    bool RequiresProcessing() const {
        return smoothing_type != "none" || taper_enabled;
    }
};

/**
 * @brief Options for state-space radiation approximation.
 * 
 * Only applies when RadiationMethod == kStateSpace.
 * 
 * The state-space method approximates the radiation impulse response as a sum
 * of exponential decay modes: K(τ) ≈ Σ H_m * exp(-α_m * τ). This enables O(1)
 * per-timestep computation vs O(N) for direct convolution.
 * 
 * These options control the fitting process that determines the number of
 * exponential modes and their parameters (H_m, α_m) from the RIRF data.
 * 
 * TUNING GUIDE:
 *   - For faster fitting: reduce max_hankel_size (e.g., 100)
 *   - For better accuracy: increase max_hankel_size (e.g., 500) and r2_threshold (e.g., 0.99)
 *   - For faster simulations: increase max_order (more modes = better fit = lower residual)
 */
struct StateSpaceOptions {
    /**
     * @brief Maximum number of state-space modes to use per DOF pair.
     * 
     * Higher values improve fit quality but increase per-step computation.
     * The fitter may use fewer modes if a good fit is achieved earlier.
     * Typical values: 4-15. Default: 10.
     */
    int max_order = 10;
    
    /**
     * @brief Minimum R² (coefficient of determination) threshold for fit acceptance.
     * 
     * The fitter tries all orders up to max_order and keeps the best R².
     * Higher orders can add fast-decaying poles that improve transient accuracy.
     * Values closer to 1.0 require better fits but may need more modes.
     * Typical values: 0.90-0.99. Default: 0.95.
     */
    double r2_threshold = 0.95;
    
    /**
     * @brief Maximum size of Hankel matrix for SVD decomposition.
     * 
     * CRITICAL PERFORMANCE PARAMETER: SVD is O(n³), so this directly controls
     * fitting time. Larger values give better accuracy but slower fitting.
     * 
     * Guidelines:
     *   - 100: Very fast fitting (~50ms), may reduce accuracy for complex kernels
     *   - 200: Good balance (default), ~300ms fitting, suitable for most cases
     *   - 500: High accuracy, ~2-5s fitting, for demanding applications
     *   - 1000: Maximum accuracy, ~20-60s fitting, rarely needed
     * 
     * If RIRF has fewer samples than this value, all samples are used.
     * Default: 200.
     */
    int max_hankel_size = 200;
    
    /**
     * @brief Number of subsamples used for R² quality check during fitting.
     * 
     * Lower values speed up fitting; higher values give more accurate R² estimates.
     * Typical values: 30-100. Default: 50.
     */
    int r2_num_samples = 50;
};

// ─────────────────────────────────────────────────────────────────────────────
// Kernel Fit Diagnostics (for state-space radiation)
// ─────────────────────────────────────────────────────────────────────────────

struct DofPairInfo {
    int dof_i = 0;
    int dof_j = 0;
    double r_squared = 0.0;
    int state_order = 0;
    int num_exp_modes = 0;
    int num_osc_modes = 0;
};

}  // namespace hydrochrono::hydro

// KernelFitDiagnostics is defined outside the namespace to avoid Eigen
// dependency in radiation_types.h when only forward-declared.
#include <Eigen/Dense>
#include <vector>

namespace hydrochrono::hydro {

struct KernelFitDiagnostics {
    std::string body_name;
    int body_index = 0;
    int num_dofs = 0;
    double min_r2 = 0.0;
    double max_r2 = 0.0;
    double mean_r2 = 0.0;
    int total_modes = 0;
    Eigen::VectorXd time;
    Eigen::MatrixXd K_actual;
    Eigen::MatrixXd K_fit;
    std::vector<DofPairInfo> dof_pairs;
};

}  // namespace hydrochrono::hydro

#endif  // HYDROC_RADIATION_TYPES_H

