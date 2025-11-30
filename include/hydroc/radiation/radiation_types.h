/*********************************************************************
 * @file  radiation_types.h
 * @brief Public types for radiation damping configuration.
 *
 * MAIN TYPES:
 *   - RadiationMethod: Selects the radiation damping calculation method
 *   - RadiationConvolutionMode: Baseline vs TaperedDirect RIRF processing
 *   - TaperedDirectOptions: Smoothing, tapering, and export settings
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
 * When using kRirfConvolution, additional options are available via
 * RadiationConvolutionMode (Baseline vs TaperedDirect).
 */
enum class RadiationMethod {
    kRirfConvolution,  ///< Direct RIRF convolution (default, existing method)
    kStateSpace        ///< State-space approximation (exponential decay modes)
};

/**
 * @brief Convolution mode for radiation damping.
 * 
 * Only applies when RadiationMethod == kRirfConvolution.
 */
enum class RadiationConvolutionMode {
    Baseline,
    TaperedDirect
};

/**
 * @brief Options for TaperedDirect RIRF preprocessing.
 * 
 * Only applies when RadiationMethod == kRirfConvolution and
 * RadiationConvolutionMode == TaperedDirect.
 */
struct TaperedDirectOptions {
    // smoothing: "sg" (Savitzky–Golay) or "moving_average"
    std::string smoothing = "sg";
    int window_length = 5;                 // odd, >= 3
    
    // RIRF truncation
    double rirf_end_time = -1.0;           // end RIRF at this time (seconds), -1.0 = use full length
    
    // Simple taper control - sensible defaults for improved stability
    double taper_start_percent = 0.8;      // start taper at 80% (taper last 20%)
    double taper_end_percent = 1.0;        // end taper at 100% of total time series  
    double taper_final_amplitude = 0.0;    // final amplitude as fraction of original (0.0 = zero, 1.0 = no change)
    bool export_plot_csv = false;          // dump before/after CSV summaries (false by default)
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
 */
struct StateSpaceOptions {
    /**
     * @brief Maximum number of exponential modes to use per DOF pair.
     * 
     * Higher values improve fit quality but increase computation time.
     * The fitter may use fewer modes if a good fit is achieved earlier.
     * Typical values: 2-10. Default: 10.
     */
    int max_order = 10;
    
    /**
     * @brief Minimum R² (coefficient of determination) threshold for fit acceptance.
     * 
     * The fitter adds modes until either R² >= r2_threshold or max_order is reached.
     * Values closer to 1.0 require better fits. Default: 0.95.
     */
    double r2_threshold = 0.95;
};

}  // namespace hydrochrono::hydro

#endif  // HYDROC_RADIATION_TYPES_H

