/*********************************************************************
 * @file  radiation_ss_model.cpp
 * @brief Implementation of RadiationStateSpaceModel.
 *********************************************************************/

#include "radiation_ss_model.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace hydrochrono::hydro {

RadiationStateSpaceModel::RadiationStateSpaceModel(int num_dofs, int num_modes)
    : num_dofs_(num_dofs), num_modes_(num_modes) {
    
    if (num_dofs <= 0) {
        throw std::invalid_argument(
            "RadiationStateSpaceModel: num_dofs must be > 0 (got " + 
            std::to_string(num_dofs) + ")");
    }
    if (num_modes <= 0) {
        throw std::invalid_argument(
            "RadiationStateSpaceModel: num_modes must be > 0 (got " + 
            std::to_string(num_modes) + ")");
    }

    // Initialize storage
    alphas_.setZero(num_modes_);
    H_.setZero(num_dofs_, num_modes_);
    z_.setZero(num_dofs_, num_modes_);
}

void RadiationStateSpaceModel::SetModeParameters(
    int mode_index, 
    double alpha, 
    const Eigen::VectorXd& h_column) {
    
    if (mode_index < 0 || mode_index >= num_modes_) {
        throw std::out_of_range(
            "RadiationStateSpaceModel::SetModeParameters: mode_index " + 
            std::to_string(mode_index) + " out of range [0, " + 
            std::to_string(num_modes_) + ")");
    }
    if (alpha <= 0.0) {
        throw std::invalid_argument(
            "RadiationStateSpaceModel::SetModeParameters: alpha must be > 0 (got " + 
            std::to_string(alpha) + ")");
    }
    if (h_column.size() != num_dofs_) {
        throw std::invalid_argument(
            "RadiationStateSpaceModel::SetModeParameters: h_column size mismatch "
            "(expected " + std::to_string(num_dofs_) + ", got " + 
            std::to_string(h_column.size()) + ")");
    }

    alphas_(mode_index) = alpha;
    H_.col(mode_index) = h_column;
}

void RadiationStateSpaceModel::Reset() {
    z_.setZero();
}

void RadiationStateSpaceModel::Step(double dt, const Eigen::VectorXd& v) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "RadiationStateSpaceModel::Step: dt must be > 0 (got " + 
            std::to_string(dt) + ")");
    }
    if (v.size() != num_dofs_) {
        throw std::invalid_argument(
            "RadiationStateSpaceModel::Step: velocity vector size mismatch "
            "(expected " + std::to_string(num_dofs_) + ", got " + 
            std::to_string(v.size()) + ")");
    }

    // Exact exponential integration for each mode:
    //   z_new = exp(-α dt) * z_old + [1 - exp(-α dt)] / α * v
    //
    // This is unconditionally stable for any α > 0 and dt > 0.
    // The coefficient [1 - exp(-α dt)] / α approaches dt as α → 0,
    // but we assume α > 0 (validated in SetModeParameters).

    for (int m = 0; m < num_modes_; ++m) {
        const double alpha = alphas_(m);
        
        // Skip uninitialized modes (alpha == 0 shouldn't happen after proper setup)
        if (alpha <= 0.0) {
            continue;
        }

        const double exp_factor = std::exp(-alpha * dt);
        const double input_coeff = (1.0 - exp_factor) / alpha;

        // Update state for mode m:
        //   z_.col(m) = exp_factor * z_.col(m) + input_coeff * v
        z_.col(m) = exp_factor * z_.col(m) + input_coeff * v;
    }
}

Eigen::VectorXd RadiationStateSpaceModel::GetForces() const {
    // Radiation force: f_rad = Σₘ Hₘ ⊙ zₘ
    //
    // With H_ and z_ both [num_dofs × num_modes]:
    //   f_rad[i] = Σₘ H_[i,m] * z_[i,m]
    //
    // This is the row-wise sum of element-wise product.

    // Element-wise product, then sum across columns (modes)
    return (H_.array() * z_.array()).rowwise().sum();
}

}  // namespace hydrochrono::hydro

