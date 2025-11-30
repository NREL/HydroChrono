/*********************************************************************
 * @file  radiation_ss_fitter.cpp
 * @brief Implementation of RadiationStateSpaceFitter.
 *
 * ALGORITHM DETAILS (matches WEC-Sim/BEMIO radiationIRFSS):
 *
 * This implementation uses a Hankel-SVD approach to identify state-space
 * parameters from an impulse response, following the WEC-Sim algorithm.
 *
 * 1. SCALE AND BUILD HANKEL MATRIX
 *    y = dt * K
 *    H = hankel(y[1:])  (Hankel matrix from scaled samples)
 *
 * 2. SVD AND ORDER SELECTION
 *    [U, S, V] = svd(H)
 *    Truncate to order O based on R² threshold
 *
 * 3. DISCRETE-TIME REALIZATION (balanced form)
 *    u1 = U[0:end-1, 0:O]
 *    u2 = U[1:end, 0:O]
 *    v1 = V[0:end-1, 0:O]
 *    sqs = sqrt(S[0:O])
 *    
 *    ubar = u1' * u2
 *    a = ubar .* ((1/sqs) * sqs')   (element-wise)
 *    b = v1[0,:] .* sqs
 *    c = u1[0,:] .* sqs'
 *    d = y[0]
 *
 * 4. BILINEAR (TUSTIN) TRANSFORM TO CONTINUOUS
 *    iidd = inv((dt/2) * (I + a))
 *    Ac = (a - I) * iidd
 *    Bc = dt * iidd * b
 *    Cc = c * iidd
 *    Dc = d - (dt/2) * c * iidd * b
 *
 * 5. KERNEL RECONSTRUCTION
 *    K_fit[k] = Cc * expm(Ac * dt * k) * Bc
 *
 *********************************************************************/

#include "radiation_ss_fitter.h"

#include <unsupported/Eigen/MatrixFunctions>
#include <algorithm>
#include <cmath>

namespace hydrochrono::hydro {

RadiationStateSpaceFitter::RadiationStateSpaceFitter(const StateSpaceOptions& opts)
    : opts_(opts) {
    // Validate options
    if (opts_.max_order <= 0) {
        opts_.max_order = 10;  // Sensible default
    }
    if (opts_.r2_threshold <= 0.0 || opts_.r2_threshold > 1.0) {
        opts_.r2_threshold = 0.95;  // Sensible default
    }
}

StateSpaceFitResult RadiationStateSpaceFitter::FitKernel(
    const Eigen::VectorXd& K,
    double dt) const {
    
    StateSpaceFitResult result;
    result.order = 0;
    result.r2 = 0.0;

    const int N = static_cast<int>(K.size());

    // Need at least 4 samples for meaningful fit
    if (N < 4) {
        return result;
    }

    // Check for trivial (zero) kernel
    const double K_norm = K.norm();
    if (K_norm < 1e-14) {
        return result;
    }

    // Validate time step
    if (dt <= 0.0) {
        return result;
    }

    // =========================================================================
    // OPTIMIZATION: Limit Hankel matrix size for faster SVD
    // The Hankel matrix only needs to capture system dynamics, not all samples.
    // User can tune via opts_.max_hankel_size (default 200).
    // =========================================================================
    const int hankel_size = std::min(N - 1, opts_.max_hankel_size);
    
    // Scale kernel
    Eigen::VectorXd y = dt * K;

    // Build reduced-size Hankel matrix
    Eigen::MatrixXd H(hankel_size, hankel_size);
    for (int i = 0; i < hankel_size; ++i) {
        for (int j = 0; j < hankel_size; ++j) {
            int idx = i + j + 1;
            H(i, j) = (idx < N) ? y(idx) : 0.0;
        }
    }

    // =========================================================================
    // OPTIMIZATION: Single SVD, compute only needed singular vectors
    // =========================================================================
    const int max_possible_order = std::min(opts_.max_order, hankel_size / 2);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::MatrixXd& U = svd.matrixU();
    const Eigen::MatrixXd& V = svd.matrixV();
    const Eigen::VectorXd& svh = svd.singularValues();

    // =========================================================================
    // Try increasing orders, reusing the same SVD decomposition
    // =========================================================================
    for (int order = 2; order <= max_possible_order; ++order) {
        StateSpaceFitResult candidate = FitFromSVD(K, dt, order, y, U, V, svh, hankel_size);
        
        if (candidate.IsValid()) {
            result = candidate;
            
            if (result.r2 >= opts_.r2_threshold) {
                break;
            }
        }
    }

    return result;
}

// Optimized method - reuses pre-computed SVD
StateSpaceFitResult RadiationStateSpaceFitter::FitFromSVD(
    const Eigen::VectorXd& K,
    double dt,
    int order,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXd& V,
    const Eigen::VectorXd& svh,
    int hankel_size) const {
    
    StateSpaceFitResult result;
    result.order = 0;
    result.r2 = 0.0;

    const int N = static_cast<int>(K.size());
    const int O = order;

    if (svh.size() < O || O < 2) {
        return result;
    }

    // Check that singular values are significant
    const double sv_threshold = 1e-10 * svh(0);
    if (svh(O - 1) < sv_threshold) {
        return result;
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Form discrete-time state-space realization (WEC-Sim style)
    // ─────────────────────────────────────────────────────────────────────────
    const int len_k = hankel_size;

    // u1 = U[0:end-1, 0:O], u2 = U[1:end, 0:O]
    Eigen::MatrixXd u1 = U.block(0, 0, len_k - 1, O);
    Eigen::MatrixXd u2 = U.block(1, 0, len_k - 1, O);
    Eigen::MatrixXd v1 = V.block(0, 0, len_k - 1, O);

    // sqs = sqrt(svh[0:O])
    Eigen::VectorXd sqs = svh.head(O).array().sqrt();

    // ubar = u1' * u2
    Eigen::MatrixXd ubar = u1.transpose() * u2;

    // a = ubar .* ((1/sqs) * sqs')
    Eigen::MatrixXd a(O, O);
    for (int i = 0; i < O; ++i) {
        for (int j = 0; j < O; ++j) {
            a(i, j) = ubar(i, j) * sqs(j) / sqs(i);
        }
    }

    // b = v1[0,:] .* sqs
    Eigen::VectorXd b = v1.row(0).transpose().array() * sqs.array();

    // c = u1[0,:] .* sqs'
    Eigen::RowVectorXd c = u1.row(0).array() * sqs.transpose().array();

    // d = y[0]
    double d = y(0);

    // ─────────────────────────────────────────────────────────────────────────
    // Bilinear (Tustin) transform: discrete to continuous
    // ─────────────────────────────────────────────────────────────────────────
    Eigen::MatrixXd I_O = Eigen::MatrixXd::Identity(O, O);
    Eigen::MatrixXd iidd = ((dt / 2.0) * (I_O + a)).inverse();

    Eigen::MatrixXd Ac = (a - I_O) * iidd;
    Eigen::VectorXd Bc = dt * (iidd * b);
    Eigen::RowVectorXd Cc = c * iidd;
    double Dc = d - (dt / 2.0) * (c * iidd * b)(0, 0);

    result.A = Ac;
    result.B = Bc;
    result.C = Cc;
    result.D = Dc;
    result.order = O;

    // ─────────────────────────────────────────────────────────────────────────
    // OPTIMIZATION: Subsample for R² computation
    // User can tune via StateSpaceOptions::r2_num_samples (default 50)
    // ─────────────────────────────────────────────────────────────────────────
    const int r2_samples = std::min(N, opts_.r2_num_samples);
    const int step = std::max(1, N / r2_samples);
    
    Eigen::VectorXd K_sub(r2_samples);
    Eigen::VectorXd K_fit_sub(r2_samples);
    
    for (int i = 0; i < r2_samples; ++i) {
        int idx = i * step;
        if (idx >= N) idx = N - 1;
        K_sub(i) = K(idx);
        
        double t = idx * dt;
        Eigen::MatrixXd expAt = (Ac * t).exp();
        K_fit_sub(i) = Cc * expAt * Bc;
    }
    
    result.r2 = ComputeR2(K_sub, K_fit_sub);

    return result;
}

Eigen::VectorXd RadiationStateSpaceFitter::ReconstructKernel(
    const StateSpaceFitResult& result,
    double dt,
    int num_samples) {
    
    if (!result.IsValid() || num_samples <= 0) {
        return Eigen::VectorXd::Zero(num_samples);
    }

    const int O = result.order;
    Eigen::VectorXd K(num_samples);

    // K[k] = C * expm(A * dt * k) * B
    // Note: For k=0, expm(0) = I, so K[0] = C * B
    for (int k = 0; k < num_samples; ++k) {
        double t = k * dt;
        Eigen::MatrixXd expAt = (result.A * t).exp();  // Matrix exponential
        K(k) = result.C * expAt * result.B;
    }

    return K;
}

double RadiationStateSpaceFitter::ComputeR2(
    const Eigen::VectorXd& K_actual,
    const Eigen::VectorXd& K_fit) {
    
    if (K_actual.size() != K_fit.size() || K_actual.size() == 0) {
        return 0.0;
    }

    // SS_res = Σ(K_actual - K_fit)²
    const double ss_res = (K_actual - K_fit).squaredNorm();

    // SS_tot = Σ(K_actual - mean(K_actual))²
    const double mean = K_actual.mean();
    const double ss_tot = (K_actual.array() - mean).matrix().squaredNorm();

    // R² = 1 - SS_res / SS_tot
    // Handle edge case where SS_tot ≈ 0 (constant kernel)
    if (ss_tot < 1e-20) {
        return (ss_res < 1e-20) ? 1.0 : 0.0;
    }

    const double r2 = 1.0 - ss_res / ss_tot;

    // Clamp to [0, 1] for numerical robustness
    return std::max(0.0, std::min(1.0, r2));
}

}  // namespace hydrochrono::hydro
