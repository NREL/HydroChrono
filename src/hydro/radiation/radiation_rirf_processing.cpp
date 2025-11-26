/*********************************************************************
 * @file  radiation_rirf_processing.cpp
 * @brief RIRF preprocessing implementation.
 *********************************************************************/

#include "radiation_rirf_processing.h"

#include <hydroc/logging.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

namespace hydrochrono::hydro {

std::vector<Eigen::Tensor<double, 3>> ProcessRirfKernels(
    int num_bodies,
    int rirf_steps,
    const Eigen::VectorXd& rirf_time_vector,
    std::function<double(int, int, int, int)> get_rirf_val,
    const TaperedDirectOptions& options,
    const std::string& diagnostics_output_dir) {
    
    const int kDofPerBody = 6;
    const int cols = kDofPerBody * num_bodies;
    const int rows = kDofPerBody;

    std::vector<Eigen::Tensor<double, 3>> processed_kernels;
    processed_kernels.resize(num_bodies);

    // SG smoothing coefficients for 5-point quadratic: [-3, 12, 17, 12, -3] / 35
    const double sg5[5] = { -3.0 / 35.0, 12.0 / 35.0, 17.0 / 35.0, 12.0 / 35.0, -3.0 / 35.0 };

    for (int b = 0; b < num_bodies; ++b) {
        Eigen::Tensor<double, 3> processed(rows, cols, rirf_steps);

        // Calculate effective steps for this body (same for all channels)
        int effective_steps = rirf_steps;
        if (options.rirf_end_time > 0.0) {
            // Calculate the step index corresponding to the end time
            double dt = rirf_time_vector[1] - rirf_time_vector[0]; // assume uniform time step
            int end_step = static_cast<int>(std::floor(options.rirf_end_time / dt));
            effective_steps = std::min(end_step, rirf_steps);
        }

        // Diagnostic aggregation across channels
        int max_tc_index = 0;
        int max_taper_len = 0;
        int max_effective_len = rirf_steps;

        for (int row_dof = 0; row_dof < rows; ++row_dof) {
            for (int col = 0; col < cols; ++col) {
                // Load raw (rho-scaled) kernel time series
                std::vector<double> k_raw(rirf_steps);
                for (int s = 0; s < rirf_steps; ++s) {
                    k_raw[s] = get_rirf_val(b, row_dof, col, s);
                }
                
                // Apply RIRF truncation if specified
                if (options.rirf_end_time > 0.0) {
                    // Truncate the raw kernel
                    k_raw.resize(effective_steps);
                }

                // Light smoothing
                std::vector<double> k_smooth(effective_steps);
                if (options.smoothing == "moving_average") {
                    const int w = std::max(3, options.window_length);
                    const int half = w / 2;
                    for (int s = 0; s < effective_steps; ++s) {
                        int a = std::max(0, s - half);
                        int b_idx = std::min(effective_steps - 1, s + half);
                        double sum = 0.0; 
                        int cnt = 0;
                        for (int i = a; i <= b_idx; ++i) { 
                            sum += k_raw[i]; 
                            ++cnt; 
                        }
                        k_smooth[s] = (cnt > 0) ? (sum / cnt) : k_raw[s];
                    }
                } else {
                    // default: SG quadratic with window 5
                    if (effective_steps >= 5) {
                        // copy edges as-is for simplicity
                        k_smooth[0] = k_raw[0];
                        k_smooth[1] = k_raw[1];
                        for (int s = 2; s <= effective_steps - 3; ++s) {
                            k_smooth[s] = sg5[0] * k_raw[s - 2] + sg5[1] * k_raw[s - 1] + sg5[2] * k_raw[s] + sg5[3] * k_raw[s + 1] + sg5[4] * k_raw[s + 2];
                        }
                        k_smooth[effective_steps - 2] = k_raw[effective_steps - 2];
                        k_smooth[effective_steps - 1] = k_raw[effective_steps - 1];
                    } else {
                        k_smooth = k_raw; // too short to smooth
                    }
                }

                // Simple taper control: start and end percentages
                int tc_index = static_cast<int>(std::floor(options.taper_start_percent * static_cast<double>(effective_steps)));
                int tc_end = static_cast<int>(std::floor(options.taper_end_percent * static_cast<double>(effective_steps)));
                tc_index = std::max(0, std::min(tc_index, effective_steps));
                tc_end = std::max(tc_index, std::min(tc_end, effective_steps));

                // Apply half-cosine taper from start to end percentage
                int taper_len = tc_end - tc_index;
                int effective_len = tc_end;

                const double pi_const = 3.14159265358979323846;
                for (int s = 0; s < effective_steps; ++s) {
                    double val = k_smooth[s];
                    if (s < tc_index) {
                        // keep original value
                    } else if (s < tc_end && taper_len > 0) {
                        double t = (static_cast<double>(s - tc_index)) / static_cast<double>(taper_len);
                        // Half-cosine taper: goes from 1.0 to taper_final_amplitude
                        // taper_final_amplitude = 0.0 means complete zero, 1.0 means no change
                        double w = options.taper_final_amplitude + (1.0 - options.taper_final_amplitude) * 0.5 * (1.0 + std::cos(pi_const * t));
                        val *= w;
                    } else {
                        val = 0.0;
                    }
                    processed(row_dof, col, s) = val;
                }
                
                // Zero out any remaining steps beyond effective_steps
                for (int s = effective_steps; s < rirf_steps; ++s) {
                    processed(row_dof, col, s) = 0.0;
                }

                if (tc_index > max_tc_index) max_tc_index = tc_index;
                if (taper_len > max_taper_len) max_taper_len = taper_len;
                if (effective_len > max_effective_len) max_effective_len = effective_len;
            }
        }

        processed_kernels[b] = std::move(processed);

        LOG_DEBUG("TaperedDirect kernel (body " << b
                  << ") Start: " << options.taper_start_percent
                  << ", End: " << options.taper_end_percent
                  << ", Final Amp: " << options.taper_final_amplitude
                  << ", RIRF End Time: " << (options.rirf_end_time > 0.0 ? std::to_string(options.rirf_end_time) + "s" : "full")
                  << ", Max Tc index: " << max_tc_index
                  << ", Max taper length: " << max_taper_len
                  << ", Max effective length: " << max_effective_len
                  << "/" << rirf_steps);

        if (options.export_plot_csv) {
            // Write small CSV summaries for inspection: times, representative channel before/after
            try {
                const std::string base = std::string("rirf_body") + std::to_string(b) + std::string("_summary.csv");
                std::filesystem::path out_dir = diagnostics_output_dir.empty() ? std::filesystem::current_path() : std::filesystem::path(diagnostics_output_dir);
                std::filesystem::path out_path = out_dir / base;
                std::ofstream ofs(out_path.string());
                ofs << "step,time,k_before,k_after\n";
                // pick row=0,col=0 as representative
                for (int s = 0; s < effective_steps; ++s) {
                    double t = (s < rirf_time_vector.size()) ? rirf_time_vector[s] : static_cast<double>(s);
                    double before = get_rirf_val(b, 0, 0, s);
                    double after = processed_kernels[b](0, 0, s);
                    ofs << s << "," << t << "," << before << "," << after << "\n";
                }
                // Only log the directory once (body 0)
                if (b == 0) {
                    LOG_INFO("RIRF CSVs written in " << out_dir.string());
                }
            } catch (...) {
                // ignore export errors
            }
        }
    }

    return processed_kernels;
}

}  // namespace hydrochrono::hydro

