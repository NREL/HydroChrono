/*********************************************************************
 * @file  radiation_ss_component.cpp
 * @brief Implementation of radiation state-space force component.
 *
 * IMPLEMENTATION NOTES:
 *
 * 1. DOF INDEXING
 *    For N bodies, we have 6N total DOFs. DOF indexing is:
 *      Body 0: DOFs 0-5 (surge, sway, heave, roll, pitch, yaw)
 *      Body 1: DOFs 6-11
 *      ...etc
 *
 * 2. RIRF INDEXING
 *    HydroData provides RIRF via GetRIRFVal(body, row_dof, col, step):
 *      - body: which body (0 to num_bodies-1)
 *      - row_dof: output DOF within body (0-5)
 *      - col: column DOF (0 to 6*num_bodies-1)
 *      - step: time step (0 to rirf_steps-1)
 *
 *    We convert to global (i,j) indexing:
 *      i = body * 6 + row_dof  (global output DOF)
 *      j = col                  (global input DOF)
 *
 * 3. FIT QUALITY
 *    We skip fitting for kernels that are essentially zero (norm < 1e-10).
 *    This avoids creating models for DOF pairs with no coupling.
 *
 * 4. RUNTIME COMPUTATION
 *    For each timestep:
 *      1. Extract velocity vector v from SystemState
 *      2. Compute dt = time - prev_time
 *      3. For each DOF pair (i,j):
 *           model.Step(dt, v[j])
 *           f_rad[i] += model.GetForce()
 *      4. Apply forces to BodyForces (with negative sign for damping)
 *********************************************************************/

#include "radiation_ss_component.h"
#include <hydroc/io/h5_reader.h>
#include <hydroc/logging.h>

#include <Eigen/Dense>
#include <algorithm>
#include <stdexcept>
#include <sstream>

namespace hydrochrono::hydro {

RadiationStateSpaceComponent::RadiationStateSpaceComponent(
    const HydroData& file_info,
    int num_bodies,
    const StateSpaceOptions& options)
    : num_bodies_(num_bodies),
      num_dofs_(kDofPerBody * num_bodies),
      options_(options),
      prev_time_(0.0),
      has_prev_time_(false) {
    
    if (num_bodies_ <= 0) {
        throw std::invalid_argument(
            "RadiationStateSpaceComponent: num_bodies must be > 0");
    }

    // Initialize state-space models from RIRF data
    InitializeFromRirf(file_info);
}

void RadiationStateSpaceComponent::InitializeFromRirf(const HydroData& file_info) {
    // Get RIRF dimensions and time vector
    const int rirf_steps = file_info.GetRIRFDims(2);
    const Eigen::VectorXd rirf_time = file_info.GetRIRFTimeVector();
    
    if (rirf_steps <= 0 || rirf_time.size() != rirf_steps) {
        throw std::runtime_error(
            "RadiationStateSpaceComponent: Invalid RIRF data dimensions");
    }

    // Compute dt from time vector (assume uniform spacing)
    const double dt = (rirf_steps > 1) 
        ? (rirf_time(rirf_time.size() - 1) - rirf_time(0)) / (rirf_steps - 1)
        : 0.01;

    // Create fitter with user options
    RadiationStateSpaceFitter fitter(options_);

    // Reserve space for models (worst case: all pairs are non-zero)
    dof_pair_models_.reserve(num_dofs_ * num_dofs_);

    int fitted_count = 0;
    int skipped_count = 0;

    // Iterate over all DOF pairs (i,j)
    for (int i = 0; i < num_dofs_; ++i) {
        for (int j = 0; j < num_dofs_; ++j) {
            // Extract kernel K_ij(t) from HydroData
            // Convert global i to (body, row_dof)
            const int body_i = i / kDofPerBody;
            const int row_dof = i % kDofPerBody;
            const int col = j;  // Column is already global index

            // Build kernel vector
            Eigen::VectorXd K(rirf_steps);
            for (int step = 0; step < rirf_steps; ++step) {
                K(step) = file_info.GetRIRFVal(body_i, row_dof, col, step);
            }

            // Skip if kernel is essentially zero (no coupling)
            const double K_norm = K.norm();
            if (K_norm < 1e-10) {
                ++skipped_count;
                continue;
            }

            // Fit the kernel
            StateSpaceFitResult result = fitter.FitKernel(K, dt);

            // Skip if fit failed
            if (!result.IsValid()) {
                ++skipped_count;
                continue;
            }

            // Create model from fit result
            RadiationStateSpaceModel model = RadiationStateSpaceModel::FromFitResult(result);
            model.Reset();

            // Store the model
            DofPairModel dof_model;
            dof_model.i = i;
            dof_model.j = j;
            dof_model.r2 = result.r2;
            dof_model.model = std::move(model);
            dof_pair_models_.push_back(std::move(dof_model));

            ++fitted_count;
        }
    }

    // Log summary using hydroc logging system
    {
        std::ostringstream oss;
        oss << "[RadiationStateSpaceComponent] Initialized: "
            << num_bodies_ << " bodies (" << num_dofs_ << " DOFs), "
            << fitted_count << " pairs fitted, "
            << skipped_count << " skipped, "
            << total_modes() << " total modes";
        if (fitted_count > 0) {
            oss << ", R² range: [" << min_r2() << ", " << max_r2() << "]";
        }
        LOG_INFO(oss.str());
    }
}

void RadiationStateSpaceComponent::Compute(
    const SystemState& state,
    double time,
    BodyForces& inout_forces) {
    
    // Validate state and forces sizes
    if (static_cast<int>(state.bodies.size()) != num_bodies_) {
        throw std::runtime_error(
            "RadiationStateSpaceComponent::Compute: state.bodies.size() mismatch");
    }
    if (static_cast<int>(inout_forces.size()) != num_bodies_) {
        throw std::runtime_error(
            "RadiationStateSpaceComponent::Compute: inout_forces.size() mismatch");
    }

    // Extract velocity vector from state
    Eigen::VectorXd v = ExtractVelocityVector(state);

    // Compute time step
    double dt;
    if (!has_prev_time_) {
        // First call - use a small default dt
        dt = 0.001;  // Will be overwritten next step
        has_prev_time_ = true;
    } else {
        dt = time - prev_time_;
    }
    prev_time_ = time;

    // Handle zero or negative dt (can happen at initialization)
    if (dt <= 0.0) {
        dt = 0.001;
    }

    // Accumulate radiation forces
    Eigen::VectorXd f_rad = Eigen::VectorXd::Zero(num_dofs_);

    for (auto& dof_model : dof_pair_models_) {
        // Step the model with velocity from input DOF j
        dof_model.model.Step(dt, v(dof_model.j));
        
        // Accumulate force contribution to output DOF i
        f_rad(dof_model.i) += dof_model.model.GetForce();
    }

    // Add forces to output (with negative sign for damping)
    AccumulateForces(f_rad, inout_forces);
}

void RadiationStateSpaceComponent::Reset() {
    for (auto& dof_model : dof_pair_models_) {
        dof_model.model.Reset();
    }
    prev_time_ = 0.0;
    has_prev_time_ = false;
}

Eigen::VectorXd RadiationStateSpaceComponent::ExtractVelocityVector(
    const SystemState& state) const {
    
    Eigen::VectorXd v(num_dofs_);
    
    for (int b = 0; b < num_bodies_; ++b) {
        const auto& body_state = state.bodies[b];
        const int offset = b * kDofPerBody;
        
        // Linear velocities (surge, sway, heave)
        v(offset + 0) = body_state.linear_velocity[0];
        v(offset + 1) = body_state.linear_velocity[1];
        v(offset + 2) = body_state.linear_velocity[2];
        
        // Angular velocities (roll, pitch, yaw)
        v(offset + 3) = body_state.angular_velocity[0];
        v(offset + 4) = body_state.angular_velocity[1];
        v(offset + 5) = body_state.angular_velocity[2];
    }
    
    return v;
}

void RadiationStateSpaceComponent::AccumulateForces(
    const Eigen::VectorXd& f_rad,
    BodyForces& inout_forces) const {
    
    // Ensure forces are properly sized
    for (int b = 0; b < num_bodies_; ++b) {
        if (static_cast<int>(inout_forces[b].size()) != kDofPerBody) {
            inout_forces[b].resize(kDofPerBody);
            inout_forces[b].setZero();
        }
    }

    // Map radiation forces to per-body forces.
    // Radiation damping opposes motion, so we SUBTRACT.
    for (int b = 0; b < num_bodies_; ++b) {
        const int offset = b * kDofPerBody;
        for (int d = 0; d < kDofPerBody; ++d) {
            inout_forces[b][d] -= f_rad(offset + d);
        }
    }
}

int RadiationStateSpaceComponent::total_modes() const {
    int total = 0;
    for (const auto& dof_model : dof_pair_models_) {
        total += dof_model.model.total_modes();
    }
    return total;
}

double RadiationStateSpaceComponent::min_r2() const {
    if (dof_pair_models_.empty()) return 0.0;
    double min_val = dof_pair_models_[0].r2;
    for (const auto& dof_model : dof_pair_models_) {
        min_val = std::min(min_val, dof_model.r2);
    }
    return min_val;
}

double RadiationStateSpaceComponent::max_r2() const {
    if (dof_pair_models_.empty()) return 0.0;
    double max_val = dof_pair_models_[0].r2;
    for (const auto& dof_model : dof_pair_models_) {
        max_val = std::max(max_val, dof_model.r2);
    }
    return max_val;
}

}  // namespace hydrochrono::hydro

