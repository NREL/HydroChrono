/*********************************************************************
 * @file  radiation_types.h
 * @brief Public types for radiation damping configuration.
 *********************************************************************/

#ifndef HYDROC_RADIATION_TYPES_H
#define HYDROC_RADIATION_TYPES_H

#include <string>

namespace hydrochrono::hydro {

/**
 * @brief Convolution mode for radiation damping.
 */
enum class RadiationConvolutionMode {
    Baseline,
    TaperedDirect
};

/**
 * @brief Options for TaperedDirect RIRF preprocessing.
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

}  // namespace hydrochrono::hydro

#endif  // HYDROC_RADIATION_TYPES_H

