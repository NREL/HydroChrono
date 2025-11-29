/*********************************************************************
 * @file  radiation_rirf_convolution.cpp
 * @brief RIRF-based radiation damping convolution implementation.
 *********************************************************************/

#include "radiation_rirf_convolution.h"

#include <algorithm>
#include <stdexcept>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hydrochrono::hydro {

const int kDofPerBody = 6;

RadiationRirfConvolution::RadiationRirfConvolution(
    int num_bodies,
    int rirf_steps,
    const Eigen::VectorXd& rirf_time_vector,
    const Eigen::VectorXd& rirf_width_vector)
    : num_bodies_(num_bodies),
      rirf_steps_(rirf_steps),
      total_dofs_(kDofPerBody * num_bodies_),
      rirf_time_vector_(rirf_time_vector),
      rirf_width_vector_(rirf_width_vector) {
    // Validate configuration parameters.
    if (num_bodies_ <= 0) {
        throw std::invalid_argument(
            "RadiationRirfConvolution: num_bodies must be > 0 (got " + 
            std::to_string(num_bodies_) + ")");
    }
    if (rirf_steps_ <= 0) {
        throw std::invalid_argument(
            "RadiationRirfConvolution: rirf_steps must be > 0 (got " + 
            std::to_string(rirf_steps_) + ")");
    }
    if (rirf_time_vector_.size() == 0) {
        throw std::invalid_argument(
            "RadiationRirfConvolution: rirf_time_vector is empty");
    }
    if (rirf_width_vector_.size() == 0) {
        throw std::invalid_argument(
            "RadiationRirfConvolution: rirf_width_vector is empty");
    }

    velocity_history_.resize(num_bodies_);
}

void RadiationRirfConvolution::RecordVelocities(
    double simulation_time,
    const std::vector<std::vector<double>>& body_velocities) {
    
    // Prevent duplicate computation within same step
    if (!time_history_.empty() && simulation_time == time_history_.front()) {
        throw std::runtime_error("Tried to compute the radiation damping convolution twice within the same time step!");
    }

    // Record current time at the front (most recent first)
    time_history_.insert(time_history_.begin(), simulation_time);

    // Record current velocities per body at the front (matching time_history_ ordering)
    for (int b = 0; b < num_bodies_; ++b) {
        auto& velocity_history_body = velocity_history_[b];
        velocity_history_body.insert(velocity_history_body.begin(), body_velocities[b]);
    }

    // Prune history older than the max RIRF time span
    const int rirf_last_index = static_cast<int>(rirf_time_vector_.size()) - 1;
    const double history_min_time = simulation_time - (rirf_last_index >= 0 ? rirf_time_vector_[rirf_last_index] : 0.0);
    PruneHistory(history_min_time);
}

void RadiationRirfConvolution::PruneHistory(double history_min_time) {
    while (time_history_.size() > 1 && time_history_[time_history_.size() - 2] < history_min_time) {
        time_history_.pop_back();
        for (int b = 0; b < num_bodies_; ++b) {
            auto& velocity_history_body = velocity_history_[b];
            if (!velocity_history_body.empty()) {
                velocity_history_body.pop_back();
            }
        }
    }
}

void RadiationRirfConvolution::InterpolateVelocity6D(
    const std::vector<std::vector<double>>& velocity_history_body,
    size_t newer_index,
    double query_time,
    double older_time,
    double newer_time,
    double out_velocity[kDofPerBody]) {
    
    if (query_time == older_time) {
        const auto& older_velocity = velocity_history_body[newer_index + 1];
        for (int dof = 0; dof < kDofPerBody; ++dof) out_velocity[dof] = older_velocity[dof];
        return;
    }
    if (query_time == newer_time) {
        const auto& newer_velocity = velocity_history_body[newer_index];
        for (int dof = 0; dof < kDofPerBody; ++dof) out_velocity[dof] = newer_velocity[dof];
        return;
    }
    if (query_time > older_time && query_time < newer_time) {
        const double time_delta   = (newer_time - older_time);
        const double weight_older = (time_delta != 0.0) ? ((newer_time - query_time) / time_delta) : 0.0;
        const double weight_newer = 1.0 - weight_older;
        const auto& older_velocity = velocity_history_body[newer_index + 1];
        const auto& newer_velocity = velocity_history_body[newer_index];
        for (int dof = 0; dof < kDofPerBody; ++dof) {
            out_velocity[dof] = weight_older * older_velocity[dof] + weight_newer * newer_velocity[dof];
        }
        return;
    }
    throw std::runtime_error("Radiation convolution: interpolation error; query_time not bracketed by history.");
}

bool RadiationRirfConvolution::AdvanceToBracket(
    const std::vector<double>& time_history,
    size_t& index,
    double query_time) {
    while ((index + 1) < time_history.size() && time_history[index + 1] > query_time) {
        ++index;
    }
    return ((index + 1) < time_history.size());
}

std::vector<double> RadiationRirfConvolution::ComputeForces(
    std::function<double(int, int, int)> get_rirf_val,
    const std::vector<Eigen::Tensor<double, 3>>* processed_rirf) {
    
    std::vector<double> force_radiation_damping(total_dofs_, 0.0);

    // Nothing to convolve with if we don't yet have at least 2 time points
    if (time_history_.size() <= 1) {
        return force_radiation_damping;
    }

    const double simulation_time = time_history_.front();
    size_t history_index = 0;  // index into descending time_history_ (front is most recent)

#ifdef _OPENMP
    const int num_threads = omp_get_max_threads();
    std::vector<std::vector<double>> thread_locals(num_threads, std::vector<double>(total_dofs_, 0.0));

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& local_out = thread_locals[tid];

        size_t history_index_local = history_index;
        #pragma omp for schedule(static)
        for (int step = 0; step < rirf_steps_; ++step) {
            const double rirf_query_time = simulation_time - rirf_time_vector_[step];

            size_t time_index = history_index_local;
            if (!AdvanceToBracket(time_history_, time_index, rirf_query_time)) {
                continue;  // not enough older history to interpolate
            }
            history_index_local = time_index;

            const double newer_time = time_history_[history_index_local];
            const double older_time = time_history_[history_index_local + 1];

            // Multibody loop: process each body's velocity history
            for (int body_index = 0; body_index < num_bodies_; ++body_index) {
                const auto& velocity_history_body = velocity_history_[body_index];
                if (velocity_history_body.size() <= history_index_local) {
                    continue;  // inconsistent; should not happen in normal flow
                }

                double interpolated_velocity_dof[kDofPerBody];
                InterpolateVelocity6D(velocity_history_body, history_index_local, rirf_query_time,
                                      older_time, newer_time, interpolated_velocity_dof);

                const double step_width = rirf_width_vector_[step];
                if (step_width == 0.0) {
                    continue;
                }
                const int body_col_offset = body_index * kDofPerBody;
                // Per-DOF loop: accumulate force contribution for each DOF
                for (int dof = 0; dof < kDofPerBody; ++dof) {
                    const int col = body_col_offset + dof;
                    const double contribution_scale = interpolated_velocity_dof[dof] * step_width;
                    if (contribution_scale == 0.0) {
                        continue;
                    }
                    for (int row = 0; row < total_dofs_; ++row) {
                        double rirf_val;
                        if (processed_rirf) {
                            // Use processed RIRF tensor (TaperedDirect mode)
                            int body_idx = row / kDofPerBody;
                            int row_dof = row % kDofPerBody;
                            rirf_val = (*processed_rirf)[body_idx](row_dof, col, step);
                        } else {
                            // Use raw RIRF via function (Baseline mode)
                            rirf_val = get_rirf_val(row, col, step);
                        }
                        local_out[row] += rirf_val * contribution_scale;
                    }
                }
            }
        }
    }

    // Deterministic combine of thread-local accumulators
    for (int t = 0; t < num_threads; ++t) {
        const auto& local = thread_locals[t];
        for (int row = 0; row < total_dofs_; ++row) {
            force_radiation_damping[row] += local[row];
        }
    }
#else
    for (int step = 0; step < rirf_steps_; ++step) {
        const double rirf_query_time = simulation_time - rirf_time_vector_[step];

        size_t time_index = history_index;
        if (!AdvanceToBracket(time_history_, time_index, rirf_query_time)) {
            break;  // not enough older history to interpolate
        }
        history_index = time_index;

        const double newer_time = time_history_[history_index];
        const double older_time = time_history_[history_index + 1];

        // Multibody loop: process each body's velocity history
        for (int body_index = 0; body_index < num_bodies_; ++body_index) {
            const auto& velocity_history_body = velocity_history_[body_index];
            if (velocity_history_body.size() <= history_index) {
                continue;  // skip if inconsistent; should not happen in normal flow
            }

            double interpolated_velocity_dof[kDofPerBody];
            InterpolateVelocity6D(velocity_history_body, history_index, rirf_query_time,
                                  older_time, newer_time, interpolated_velocity_dof);

            const double step_width = rirf_width_vector_[step];
            if (step_width == 0.0) {
                continue;
            }
            const int body_col_offset = body_index * kDofPerBody;
            // Per-DOF loop: accumulate force contribution for each DOF
            for (int dof = 0; dof < kDofPerBody; ++dof) {
                const int col = body_col_offset + dof;
                const double contribution_scale = interpolated_velocity_dof[dof] * step_width;
                if (contribution_scale == 0.0) {
                    continue;
                }
                for (int row = 0; row < total_dofs_; ++row) {
                    double rirf_val;
                    if (processed_rirf) {
                        // Use processed RIRF tensor (TaperedDirect mode)
                        int body_idx = row / kDofPerBody;
                        int row_dof = row % kDofPerBody;
                        rirf_val = (*processed_rirf)[body_idx](row_dof, col, step);
                    } else {
                        // Use raw RIRF via function (Baseline mode)
                        rirf_val = get_rirf_val(row, col, step);
                    }
                    force_radiation_damping[row] += rirf_val * contribution_scale;
                }
            }
        }
    }
#endif

    return force_radiation_damping;
}

}  // namespace hydrochrono::hydro

