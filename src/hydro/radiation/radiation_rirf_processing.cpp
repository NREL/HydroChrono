/*********************************************************************
 * @file  radiation_rirf_processing.cpp
 * @brief RIRF kernel processing implementation (smoothing + tapering).
 *
 * Truncation is applied upstream by HydroSystem::CreateRadiationComponent()
 * before this module is invoked.  See radiation_rirf_processing.h for the
 * full processing pipeline description.
 *********************************************************************/

#include "radiation_rirf_processing.h"

#include <hydroc/logging.h>
#include <hydroc/math_constants.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <array>
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
    const RadiationKernelProcessing& options,
    const std::string& diagnostics_output_dir) {
    
    const int kDofPerBody = 6;
    const int cols = kDofPerBody * num_bodies;
    const int rows = kDofPerBody;

    std::vector<Eigen::Tensor<double, 3>> processed_kernels;
    processed_kernels.resize(num_bodies);

    // Savitzky-Golay 5-point quadratic smoothing coefficients: [-3, 12, 17, 12, -3] / 35.
    // Reference: Savitzky & Golay (1964), Analytical Chemistry 36(8), pp. 1627-1639.
    static constexpr std::array<double, 5> kSG5 = {
        -3.0 / 35.0, 12.0 / 35.0, 17.0 / 35.0, 12.0 / 35.0, -3.0 / 35.0
    };

    const bool do_smoothing = (options.smoothing_type != "none");
    const bool do_taper = options.taper_enabled;

    for (int b = 0; b < num_bodies; ++b) {
        Eigen::Tensor<double, 3> processed(rows, cols, rirf_steps);

        for (int row_dof = 0; row_dof < rows; ++row_dof) {
            for (int col = 0; col < cols; ++col) {
                std::vector<double> k_raw(rirf_steps);
                for (int s = 0; s < rirf_steps; ++s) {
                    k_raw[s] = get_rirf_val(b, row_dof, col, s);
                }

                // Step 1: Smoothing (if configured)
                std::vector<double> k_smooth(rirf_steps);
                if (do_smoothing && options.smoothing_type == "moving_average") {
                    const int w = std::max(3, options.smoothing_window);
                    const int half = w / 2;
                    for (int s = 0; s < rirf_steps; ++s) {
                        int a = std::max(0, s - half);
                        int b_idx = std::min(rirf_steps - 1, s + half);
                        double sum = 0.0; 
                        int cnt = 0;
                        for (int i = a; i <= b_idx; ++i) { 
                            sum += k_raw[i]; 
                            ++cnt; 
                        }
                        k_smooth[s] = (cnt > 0) ? (sum / cnt) : k_raw[s];
                    }
                } else if (do_smoothing) {
                    // Default: Savitzky-Golay quadratic, window = 5.
                    // Edge points are copied unchanged (filter needs ±2 neighbours).
                    if (rirf_steps >= 5) {
                        k_smooth[0] = k_raw[0];
                        k_smooth[1] = k_raw[1];
                        for (int s = 2; s <= rirf_steps - 3; ++s) {
                            k_smooth[s] = kSG5[0] * k_raw[s - 2]
                                        + kSG5[1] * k_raw[s - 1]
                                        + kSG5[2] * k_raw[s]
                                        + kSG5[3] * k_raw[s + 1]
                                        + kSG5[4] * k_raw[s + 2];
                        }
                        k_smooth[rirf_steps - 2] = k_raw[rirf_steps - 2];
                        k_smooth[rirf_steps - 1] = k_raw[rirf_steps - 1];
                    } else {
                        k_smooth = k_raw;
                    }
                } else {
                    k_smooth = k_raw;
                }

                // Step 2: Tapering (if configured).
                // Half-cosine window: w(t) ramps from 1 to taper_final_amplitude
                // over [taper_start_percent, taper_end_percent] of the RIRF length,
                // then zeros beyond taper_end_percent.
                if (do_taper) {
                    const int tc_start = std::clamp(
                        static_cast<int>(std::floor(options.taper_start_percent * rirf_steps)),
                        0, rirf_steps);
                    const int tc_end = std::clamp(
                        static_cast<int>(std::floor(options.taper_end_percent * rirf_steps)),
                        tc_start, rirf_steps);
                    const int taper_len = tc_end - tc_start;

                    for (int s = 0; s < rirf_steps; ++s) {
                        double val = k_smooth[s];
                        if (s >= tc_start && s < tc_end && taper_len > 0) {
                            const double frac = static_cast<double>(s - tc_start)
                                              / static_cast<double>(taper_len);
                            const double w = options.taper_final_amplitude
                                           + (1.0 - options.taper_final_amplitude)
                                             * 0.5 * (1.0 + std::cos(M_PI * frac));
                            val *= w;
                        } else if (s >= tc_end) {
                            val = 0.0;
                        }
                        processed(row_dof, col, s) = val;
                    }
                } else {
                    for (int s = 0; s < rirf_steps; ++s) {
                        processed(row_dof, col, s) = k_smooth[s];
                    }
                }
            }
        }

        processed_kernels[b] = std::move(processed);

        LOG_DEBUG("RIRF kernel processing (body " << b
                  << ") Smoothing: " << options.smoothing_type
                  << ", Taper: " << (do_taper ? "enabled" : "disabled")
                  << (do_taper ? (", Start: " + std::to_string(options.taper_start_percent)
                                + ", End: " + std::to_string(options.taper_end_percent)
                                + ", Final Amp: " + std::to_string(options.taper_final_amplitude))
                              : "")
                  << ", Steps: " << rirf_steps);

        if (options.export_csv) {
            try {
                const std::string base = std::string("rirf_body") + std::to_string(b) + std::string("_summary.csv");
                std::filesystem::path out_dir = diagnostics_output_dir.empty() ? std::filesystem::current_path() : std::filesystem::path(diagnostics_output_dir);
                std::filesystem::path out_path = out_dir / base;
                std::ofstream ofs(out_path.string());
                ofs << "step,time,k_before,k_after\n";
                for (int s = 0; s < rirf_steps; ++s) {
                    double t = (s < rirf_time_vector.size()) ? rirf_time_vector[s] : static_cast<double>(s);
                    double before = get_rirf_val(b, 0, 0, s);
                    double after = processed_kernels[b](0, 0, s);
                    ofs << s << "," << t << "," << before << "," << after << "\n";
                }
                if (b == 0) {
                    LOG_INFO("RIRF CSVs written in " << out_dir.string());
                }
            } catch (const std::exception& e) {
                LOG_WARNING("RIRF CSV export failed for body " << b << ": " << e.what());
            } catch (...) {
                LOG_WARNING("RIRF CSV export failed for body " << b << " (unknown error)");
            }
        }
    }

    return processed_kernels;
}

}  // namespace hydrochrono::hydro
