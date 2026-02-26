/*********************************************************************
 * @file  radiation_rirf_processing.h
 * @brief RIRF kernel processing: smoothing, tapering, CSV export.
 *
 * Truncation is applied before this module is called.  Specifically,
 * HydroSystem::CreateRadiationComponent() trims the RIRF time/width
 * vectors to [0, T] when SetRadiationTruncationTime(T) has been set.
 * This module then applies optional smoothing and tapering to the
 * already-truncated data.
 *********************************************************************/

#ifndef HYDRO_RADIATION_RIRF_PROCESSING_H
#define HYDRO_RADIATION_RIRF_PROCESSING_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <string>
#include <vector>

#include <hydroc/radiation/radiation_types.h>

namespace hydrochrono::hydro {

/**
 * @brief Apply smoothing and/or tapering to RIRF kernels.
 * 
 * Expects pre-truncated inputs (see file-level comment). Applies:
 *   1. Smoothing (Savitzky-Golay or moving average) if configured
 *   2. Half-cosine tapering window if configured
 * 
 * @param num_bodies Number of bodies in the system
 * @param rirf_steps Number of time steps in the (possibly truncated) RIRF
 * @param rirf_time_vector Time vector for RIRF (used for CSV export timestamps)
 * @param get_rirf_val Function to get raw RIRF value: (body, row_dof, col, step) -> value
 * @param options Kernel processing options (smoothing + tapering)
 * @param diagnostics_output_dir Directory for CSV exports (empty = current directory)
 * 
 * @return Processed RIRF tensors, one per body [dof x col x step]
 */
std::vector<Eigen::Tensor<double, 3>> ProcessRirfKernels(
    int num_bodies,
    int rirf_steps,
    const Eigen::VectorXd& rirf_time_vector,
    std::function<double(int, int, int, int)> get_rirf_val,
    const RadiationKernelProcessing& options,
    const std::string& diagnostics_output_dir = "");

}  // namespace hydrochrono::hydro

#endif  // HYDRO_RADIATION_RIRF_PROCESSING_H

