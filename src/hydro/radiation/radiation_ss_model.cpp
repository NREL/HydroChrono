/*********************************************************************
 * @file  radiation_ss_model.cpp
 * @brief Implementation of RadiationStateSpaceModel.
 *
 * EXACT INTEGRATION FORMULAS:
 *
 * 1. EXPONENTIAL MODE (z' = -α*z + b*v)
 *    With constant v over [t, t+dt]:
 *    z(t+dt) = exp(-α*dt) * z(t) + b * (1 - exp(-α*dt)) / α * v
 *
 * 2. OSCILLATORY MODE
 *    [z_c]' = [-α  -ω] [z_c] + [b_c] * v
 *    [z_s]    [ ω  -α] [z_s]   [b_s]
 *
 *    The homogeneous solution uses rotation + decay:
 *    exp(A*t) = exp(-α*t) * [cos(ω*t)  -sin(ω*t)]
 *                           [sin(ω*t)   cos(ω*t)]
 *
 *    The particular integral for constant v involves:
 *    ∫[0,dt] exp(-α*s) * cos(ω*s) ds
 *    ∫[0,dt] exp(-α*s) * sin(ω*s) ds
 *
 *********************************************************************/

#include "radiation_ss_model.h"
#include "radiation_ss_fitter.h"

#include <cmath>
#include <complex>
#include <algorithm>
#include <stdexcept>

namespace hydrochrono::hydro {

RadiationStateSpaceModel RadiationStateSpaceModel::FromFitResult(const StateSpaceFitResult& result) {
    RadiationStateSpaceModel model;

    if (!result.IsValid()) {
        return model;
    }

    const int n = result.order;
    const Eigen::MatrixXd& A = result.A;
    const Eigen::VectorXd& B = result.B;
    const Eigen::RowVectorXd& C = result.C;

    // Compute eigendecomposition of A
    Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
    const Eigen::VectorXcd& eigenvalues = solver.eigenvalues();
    const Eigen::MatrixXcd& eigenvectors = solver.eigenvectors();

    // Compute left eigenvectors (rows of V^{-1})
    Eigen::MatrixXcd V_inv = eigenvectors.inverse();

    // Track which eigenvalues we've processed (for conjugate pairs)
    std::vector<bool> processed(n, false);

    for (int k = 0; k < n; ++k) {
        if (processed[k]) continue;

        std::complex<double> lambda = eigenvalues(k);
        double real_part = lambda.real();
        double imag_part = lambda.imag();

        // Check for stability (real part should be negative for decay)
        if (real_part >= 0) {
            processed[k] = true;
            continue;  // Skip unstable modes
        }

        if (std::abs(imag_part) < 1e-10) {
            // Real eigenvalue -> pure exponential mode
            double alpha = -real_part;  // α > 0 for decay

            // Compute input/output gains
            std::complex<double> b_complex = V_inv.row(k) * B.cast<std::complex<double>>();
            std::complex<double> H_complex = C.cast<std::complex<double>>() * eigenvectors.col(k);

            double b = b_complex.real();
            double H = H_complex.real();

            model.AddExponentialMode(alpha, b, H);
            processed[k] = true;
        } else {
            // Complex eigenvalue -> look for conjugate pair
            int conj_idx = -1;
            for (int j = k + 1; j < n; ++j) {
                if (!processed[j]) {
                    std::complex<double> other = eigenvalues(j);
                    if (std::abs(other.real() - real_part) < 1e-10 &&
                        std::abs(other.imag() + imag_part) < 1e-10) {
                        conj_idx = j;
                        break;
                    }
                }
            }

            if (conj_idx < 0) {
                // No conjugate found, skip
                processed[k] = true;
                continue;
            }

            // Process conjugate pair
            double alpha = -real_part;  // α > 0
            double omega = std::abs(imag_part);  // ω > 0

            // For conjugate pair, compute the real representation
            // K(t) = 2*Re[g * exp(λ*t)] where g = (C*r_k) * (l_k*B)
            //      = 2*exp(-α*t) * [Re(g)*cos(ωt) - Im(g)*sin(ωt)]
            std::complex<double> l_k_B = V_inv.row(k) * B.cast<std::complex<double>>();
            std::complex<double> C_r_k = C.cast<std::complex<double>>() * eigenvectors.col(k);
            std::complex<double> g = C_r_k * l_k_B;

            double H_c = 2.0 * g.real();
            double H_s = -2.0 * g.imag();

            // Standard representation with b_c=1, b_s=0
            double b_c = 1.0;
            double b_s = 0.0;

            model.AddOscillatoryMode(alpha, omega, b_c, b_s, H_c, H_s);

            processed[k] = true;
            processed[conj_idx] = true;
        }
    }

    return model;
}

void RadiationStateSpaceModel::AddExponentialMode(double alpha, double b, double H) {
    if (alpha <= 0) {
        throw std::invalid_argument("alpha must be positive");
    }
    exp_modes_.emplace_back(alpha, b, H);
}

void RadiationStateSpaceModel::AddOscillatoryMode(double alpha, double omega,
                                                   double b_c, double b_s,
                                                   double H_c, double H_s) {
    if (alpha <= 0) {
        throw std::invalid_argument("alpha must be positive");
    }
    if (omega <= 0) {
        throw std::invalid_argument("omega must be positive");
    }
    OscillatoryMode mode;
    mode.alpha = alpha;
    mode.omega = omega;
    mode.b_c = b_c;
    mode.b_s = b_s;
    mode.H_c = H_c;
    mode.H_s = H_s;
    mode.z_c = 0.0;
    mode.z_s = 0.0;
    osc_modes_.push_back(mode);
}

void RadiationStateSpaceModel::Reset() {
    for (auto& mode : exp_modes_) {
        mode.z = 0.0;
    }
    for (auto& mode : osc_modes_) {
        mode.z_c = 0.0;
        mode.z_s = 0.0;
    }
}

void RadiationStateSpaceModel::Step(double dt, double v) {
    if (dt <= 0) {
        throw std::invalid_argument("dt must be positive");
    }

    // Update exponential modes
    // z(t+dt) = exp(-α*dt) * z(t) + b * (1 - exp(-α*dt)) / α * v
    for (auto& mode : exp_modes_) {
        double exp_decay = std::exp(-mode.alpha * dt);
        double gain = mode.b * (1.0 - exp_decay) / mode.alpha;
        mode.z = exp_decay * mode.z + gain * v;
    }

    // Update oscillatory modes using rotation + decay + particular integral
    for (auto& mode : osc_modes_) {
        double alpha = mode.alpha;
        double omega = mode.omega;
        double exp_decay = std::exp(-alpha * dt);
        double cos_wt = std::cos(omega * dt);
        double sin_wt = std::sin(omega * dt);

        // Homogeneous solution (rotation + decay)
        double z_c_old = mode.z_c;
        double z_s_old = mode.z_s;
        double z_c_hom = exp_decay * (cos_wt * z_c_old - sin_wt * z_s_old);
        double z_s_hom = exp_decay * (sin_wt * z_c_old + cos_wt * z_s_old);

        // Particular integral for constant input v
        double denom = alpha * alpha + omega * omega;
        double int_cos = (alpha - exp_decay * (alpha * cos_wt - omega * sin_wt)) / denom;
        double int_sin = (omega - exp_decay * (omega * cos_wt + alpha * sin_wt)) / denom;

        // Contribution from input [b_c; b_s] * v
        double z_c_part = (mode.b_c * int_cos - mode.b_s * int_sin) * v;
        double z_s_part = (mode.b_c * int_sin + mode.b_s * int_cos) * v;

        mode.z_c = z_c_hom + z_c_part;
        mode.z_s = z_s_hom + z_s_part;
    }
}

double RadiationStateSpaceModel::GetForce() const {
    double force = 0.0;

    // Contribution from exponential modes
    for (const auto& mode : exp_modes_) {
        force += mode.H * mode.z;
    }

    // Contribution from oscillatory modes
    for (const auto& mode : osc_modes_) {
        force += mode.H_c * mode.z_c + mode.H_s * mode.z_s;
    }

    return force;
}

Eigen::VectorXd RadiationStateSpaceModel::ReconstructKernel(double dt, int num_samples) const {
    Eigen::VectorXd K(num_samples);

    for (int k = 0; k < num_samples; ++k) {
        double t = k * dt;
        double val = 0.0;

        // Exponential mode contributions: K(t) = H * b * exp(-α*t)
        for (const auto& mode : exp_modes_) {
            val += mode.H * mode.b * std::exp(-mode.alpha * t);
        }

        // Oscillatory mode contributions
        // Impulse response: z_c(t), z_s(t) starting from z_c(0)=b_c, z_s(0)=b_s
        for (const auto& mode : osc_modes_) {
            double exp_decay = std::exp(-mode.alpha * t);
            double cos_wt = std::cos(mode.omega * t);
            double sin_wt = std::sin(mode.omega * t);

            double z_c_t = exp_decay * (mode.b_c * cos_wt - mode.b_s * sin_wt);
            double z_s_t = exp_decay * (mode.b_c * sin_wt + mode.b_s * cos_wt);

            val += mode.H_c * z_c_t + mode.H_s * z_s_t;
        }

        K(k) = val;
    }

    return K;
}

}  // namespace hydrochrono::hydro
