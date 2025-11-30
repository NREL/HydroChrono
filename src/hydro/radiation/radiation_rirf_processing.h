/*********************************************************************
 * @file  radiation_rirf_processing.h
 * @brief RIRF preprocessing: truncation, smoothing, tapering, CSV export.
 *********************************************************************/

#ifndef HYDRO_RADIATION_RIRF_PROCESSING_H
#define HYDRO_RADIATION_RIRF_PROCESSING_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <string>
#include <vector>

// Include canonical TaperedDirectOptions from public header
#include <hydroc/radiation/radiation_types.h>

namespace hydrochrono::hydro {

/**
 * @brief Preprocess RIRF kernels for TaperedDirect mode.
 * 
 * Applies truncation, smoothing (moving average or Savitzky-Golay), and tapering
 * to raw RIRF kernels. Optionally exports CSV diagnostics.
 * 
 * @param num_bodies Number of bodies in the system
 * @param rirf_steps Number of time steps in RIRF
 * @param rirf_time_vector Time vector for RIRF (for truncation and CSV export)
 * @param get_rirf_val Function to get raw RIRF value: (body, row_dof, col, step) -> value
 * @param options Preprocessing options
 * @param diagnostics_output_dir Directory for CSV exports (empty = current directory)
 * 
 * @return Processed RIRF tensors, one per body [dof x col x step]
 */
std::vector<Eigen::Tensor<double, 3>> ProcessRirfKernels(
    int num_bodies,
    int rirf_steps,
    const Eigen::VectorXd& rirf_time_vector,
    std::function<double(int, int, int, int)> get_rirf_val,
    const TaperedDirectOptions& options,
    const std::string& diagnostics_output_dir = "");

}  // namespace hydrochrono::hydro

#endif  // HYDRO_RADIATION_RIRF_PROCESSING_H

