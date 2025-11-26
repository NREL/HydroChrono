/*********************************************************************
 * @file  radiation_model.cpp
 * @brief Implementation of radiation damping force model.
 *********************************************************************/

#include "radiation_model.h"
#include <hydroc/h5fileinfo.h>

#include <Eigen/Dense>
#include <cassert>
#include <algorithm>
#include <stdexcept>

namespace hydrochrono::hydro {

RadiationModel::RadiationModel(
    const HydroData& file_info,
    int num_bodies,
    int rirf_steps,
    const Eigen::VectorXd& rirf_time_vector,
    const Eigen::VectorXd& rirf_width_vector,
    RadiationConvolutionMode convolution_mode,
    const TaperedDirectOptions& tapered_opts,
    const std::string& diagnostics_output_dir)
    : file_info_(file_info),
      num_bodies_(num_bodies),
      convolution_mode_(convolution_mode),
      tapered_opts_(tapered_opts),
      diagnostics_output_dir_(diagnostics_output_dir),
      convolution_(std::make_unique<RadiationRirfConvolution>(
          num_bodies, rirf_steps, rirf_time_vector, rirf_width_vector)),
      rirf_time_vector_(rirf_time_vector),
      rirf_width_vector_(rirf_width_vector),
      rirf_processed_ready_(false) {
    assert(num_bodies_ > 0);
}

void RadiationModel::EnsureProcessedRIRF() {
    if (rirf_processed_ready_) {
        return;
    }

    const int steps = file_info_.GetRIRFDims(2);
    
    // Create lambda to get RIRF values from file_info_
    auto get_rirf_val = [this](int body, int row_dof, int col, int step) -> double {
        return file_info_.GetRIRFVal(body, row_dof, col, step);
    };

    // Process RIRF kernels using the dedicated module
    rirf_processed_ = ProcessRirfKernels(
        num_bodies_, steps, rirf_time_vector_, get_rirf_val, tapered_opts_, diagnostics_output_dir_);

    rirf_processed_ready_ = true;
}

double RadiationModel::GetRIRFval(int row, int col, int st) {
    if (row < 0 || row >= kDofPerBody * num_bodies_ || col < 0 || col >= kDofPerBody * num_bodies_ || st < 0 ||
        st >= file_info_.GetRIRFDims(2)) {
        throw std::out_of_range("rirfval index out of range in RadiationModel");
    }

    int body_index = row / kDofPerBody;
    int col_dof    = col % kDofPerBody;
    int row_dof    = row % kDofPerBody;

    if (convolution_mode_ == RadiationConvolutionMode::TaperedDirect) {
        EnsureProcessedRIRF();
        // processed tensor is scaled by rho already
        const auto& tensor = rirf_processed_[body_index];
        return tensor(row_dof, col, st);
    }

    return file_info_.GetRIRFVal(body_index, row_dof, col, st);
}

void RadiationModel::Compute(const SystemState& state,
                             double time,
                             BodyForces& inout_forces) {
    assert(static_cast<int>(state.bodies.size()) == num_bodies_);
    assert(static_cast<int>(inout_forces.size()) == num_bodies_);

    const int rirf_steps = file_info_.GetRIRFDims(2);
    const int total_dofs = kDofPerBody * num_bodies_;

    assert(total_dofs > 0 && rirf_steps > 0);

    // If using TaperedDirect, ensure processed kernel is ready
    if (convolution_mode_ == RadiationConvolutionMode::TaperedDirect) {
        EnsureProcessedRIRF();
    }

    // Prepare body velocities for recording
    std::vector<std::vector<double>> body_velocities(num_bodies_);
    for (int b = 0; b < num_bodies_; ++b) {
        const auto& body_state = state.bodies[b];
        body_velocities[b] = {
            body_state.linear_velocity[0], body_state.linear_velocity[1], body_state.linear_velocity[2],
            body_state.angular_velocity[0], body_state.angular_velocity[1], body_state.angular_velocity[2]
        };
    }

    // Record velocities in the convolution module
    convolution_->RecordVelocities(time, body_velocities);

    // Create lambda to get RIRF values (for Baseline mode)
    auto get_rirf_val = [this](int row, int col, int step) -> double {
        return GetRIRFval(row, col, step);
    };

    // Compute forces using the convolution module
    const std::vector<Eigen::Tensor<double, 3>>* processed_rirf = 
        (convolution_mode_ == RadiationConvolutionMode::TaperedDirect && rirf_processed_ready_) 
        ? &rirf_processed_ : nullptr;

    std::vector<double> force_flat = convolution_->ComputeForces(get_rirf_val, processed_rirf);

    // Convert flat 6N vector to BodyForces (per-body 6-DOF vectors)
    // Ensure forces are properly sized
    for (int b = 0; b < num_bodies_; ++b) {
        if (static_cast<int>(inout_forces[b].size()) != kDofPerBody) {
            inout_forces[b].resize(kDofPerBody);
            inout_forces[b].setZero();
        }
    }

    // Map flat vector to per-body forces and add to inout_forces
    for (int b = 0; b < num_bodies_; ++b) {
        const int body_offset = kDofPerBody * b;
        for (int i = 0; i < kDofPerBody; ++i) {
            inout_forces[b][i] += force_flat[body_offset + i];
        }
    }
}

}  // namespace hydrochrono::hydro

