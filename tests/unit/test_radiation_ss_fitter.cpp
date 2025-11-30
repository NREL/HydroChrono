/*********************************************************************
 * @file  test_radiation_ss_fitter.cpp
 * @brief Unit tests for RadiationStateSpaceFitter.
 *
 * These tests verify that the fitter can recover state-space models
 * from synthetic kernels (both exponential and oscillatory).
 *
 * TEST SCENARIOS:
 *
 * 1. Single Exponential Kernel
 *    K(t) = H * exp(-α t)
 *    The fitter should achieve high R².
 *
 * 2. Two-Mode Exponential Kernel
 *    K(t) = H1 * exp(-α1 t) + H2 * exp(-α2 t)
 *    The fitter should find order ≥ 2 and good R².
 *
 * 3. Oscillatory Kernel
 *    K(t) = A * exp(-α t) * cos(ω t)
 *    The fitter should achieve high R² with complex eigenvalues.
 *
 * 4. Zero/Trivial Kernel
 *    K(t) = 0
 *    The fitter should return order = 0 gracefully.
 *
 * 5. Integration with RadiationStateSpaceModel
 *    Fit a kernel, use the result in the model, verify forces match.
 *********************************************************************/

#include "../../src/hydro/radiation/radiation_ss_fitter.h"
#include "../../src/hydro/radiation/radiation_ss_model.h"

#include <cmath>
#include <iostream>
#include <iomanip>

// Simple test framework macros
#define TEST_ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "FAILED: " << message << std::endl; \
            std::cerr << "  at " << __FILE__ << ":" << __LINE__ << std::endl; \
            return false; \
        } \
    } while (0)

#define TEST_ASSERT_NEAR(actual, expected, tolerance, message) \
    do { \
        double _actual = (actual); \
        double _expected = (expected); \
        double _diff = std::abs(_actual - _expected); \
        if (_diff > (tolerance)) { \
            std::cerr << "FAILED: " << message << std::endl; \
            std::cerr << "  Expected: " << _expected << ", Got: " << _actual \
                      << ", Diff: " << _diff << ", Tolerance: " << (tolerance) << std::endl; \
            std::cerr << "  at " << __FILE__ << ":" << __LINE__ << std::endl; \
            return false; \
        } \
    } while (0)

using namespace hydrochrono::hydro;

// =============================================================================
// Helper: Generate synthetic kernels
// =============================================================================
Eigen::VectorXd GenerateSingleExpKernel(double H, double alpha, double dt, int N) {
    Eigen::VectorXd K(N);
    for (int k = 0; k < N; ++k) {
        K(k) = H * std::exp(-alpha * k * dt);
    }
    return K;
}

Eigen::VectorXd GenerateTwoExpKernel(double H1, double alpha1, 
                                      double H2, double alpha2, 
                                      double dt, int N) {
    Eigen::VectorXd K(N);
    for (int k = 0; k < N; ++k) {
        double t = k * dt;
        K(k) = H1 * std::exp(-alpha1 * t) + H2 * std::exp(-alpha2 * t);
    }
    return K;
}

Eigen::VectorXd GenerateOscillatoryKernel(double A, double alpha, double omega, 
                                           double dt, int N) {
    Eigen::VectorXd K(N);
    for (int k = 0; k < N; ++k) {
        double t = k * dt;
        K(k) = A * std::exp(-alpha * t) * std::cos(omega * t);
    }
    return K;
}

// =============================================================================
// Test: Single exponential kernel recovery
// =============================================================================
bool TestSingleExponential() {
    std::cout << "  Testing single exponential kernel..." << std::endl;

    const double H_true = 5.0;
    const double alpha_true = 2.0;
    const double dt = 0.01;
    const int N = 500;

    Eigen::VectorXd K = GenerateSingleExpKernel(H_true, alpha_true, dt, N);

    StateSpaceOptions opts;
    opts.max_order = 5;
    opts.r2_threshold = 0.99;

    RadiationStateSpaceFitter fitter(opts);
    StateSpaceFitResult result = fitter.FitKernel(K, dt);

    TEST_ASSERT(result.IsValid(), "Fit should be valid for single exponential");
    TEST_ASSERT(result.order >= 1, "Should find at least 1 mode");
    TEST_ASSERT(result.r2 >= 0.99, "R² should be >= 0.99 for clean exponential");

    // Verify by creating model and reconstructing kernel
    RadiationStateSpaceModel model = RadiationStateSpaceModel::FromFitResult(result);
    Eigen::VectorXd K_fit = model.ReconstructKernel(dt, N);

    double max_rel_error = 0.0;
    for (int k = 0; k < N; ++k) {
        if (std::abs(K(k)) > 0.001) {
            double rel_error = std::abs(K(k) - K_fit(k)) / std::abs(K(k));
            max_rel_error = std::max(max_rel_error, rel_error);
        }
    }

    TEST_ASSERT(max_rel_error < 0.01, "Kernel reconstruction error should be < 1%");

    std::cout << "    R² = " << std::fixed << std::setprecision(6) << result.r2 << std::endl;
    std::cout << "    Order = " << result.order << std::endl;
    std::cout << "    Model: " << model.num_exp_modes() << " exp + " 
              << model.num_osc_modes() << " osc modes" << std::endl;
    std::cout << "    Max relative error = " << max_rel_error << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Two-mode exponential kernel
// =============================================================================
bool TestTwoModeExponential() {
    std::cout << "  Testing two-mode exponential kernel..." << std::endl;

    const double H1_true = 3.0, alpha1_true = 1.0;
    const double H2_true = 2.0, alpha2_true = 5.0;
    const double dt = 0.01;
    const int N = 500;

    Eigen::VectorXd K = GenerateTwoExpKernel(H1_true, alpha1_true, 
                                              H2_true, alpha2_true, 
                                              dt, N);

    StateSpaceOptions opts;
    opts.max_order = 10;
    opts.r2_threshold = 0.99;

    RadiationStateSpaceFitter fitter(opts);
    StateSpaceFitResult result = fitter.FitKernel(K, dt);

    TEST_ASSERT(result.IsValid(), "Fit should be valid for two-mode kernel");
    TEST_ASSERT(result.order >= 2, "Should find at least 2 modes");
    TEST_ASSERT(result.r2 >= 0.98, "R² should be >= 0.98");

    // Verify reconstruction
    Eigen::VectorXd K_fit = RadiationStateSpaceFitter::ReconstructKernel(result, dt, N);

    double max_rel_error = 0.0;
    for (int k = 0; k < N; ++k) {
        if (std::abs(K(k)) > 0.01) {
            double rel_error = std::abs(K(k) - K_fit(k)) / std::abs(K(k));
            max_rel_error = std::max(max_rel_error, rel_error);
        }
    }

    TEST_ASSERT(max_rel_error < 0.1, "Max relative error should be < 10%");

    std::cout << "    R² = " << std::fixed << std::setprecision(6) << result.r2 << std::endl;
    std::cout << "    Order = " << result.order << std::endl;
    std::cout << "    Max relative error = " << max_rel_error << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Oscillatory kernel (damped cosine)
// =============================================================================
bool TestOscillatoryKernel() {
    std::cout << "  Testing oscillatory kernel..." << std::endl;

    const double A_true = 5.0;
    const double alpha_true = 0.8;
    const double omega_true = 3.0;
    const double dt = 0.02;
    const int N = 500;

    Eigen::VectorXd K = GenerateOscillatoryKernel(A_true, alpha_true, omega_true, dt, N);

    StateSpaceOptions opts;
    opts.max_order = 10;
    opts.r2_threshold = 0.999;

    RadiationStateSpaceFitter fitter(opts);
    StateSpaceFitResult result = fitter.FitKernel(K, dt);

    TEST_ASSERT(result.IsValid(), "Fit should be valid for oscillatory kernel");
    TEST_ASSERT(result.r2 >= 0.999, "R² should be >= 0.999 for clean oscillatory");

    // Create model and verify
    RadiationStateSpaceModel model = RadiationStateSpaceModel::FromFitResult(result);
    TEST_ASSERT(model.num_osc_modes() >= 1, "Should have at least 1 oscillatory mode");

    Eigen::VectorXd K_fit = model.ReconstructKernel(dt, N);

    double max_error = 0.0;
    for (int k = 0; k < N; ++k) {
        double error = std::abs(K(k) - K_fit(k));
        max_error = std::max(max_error, error);
    }

    TEST_ASSERT(max_error < 0.1, "Max absolute error should be < 0.1");

    std::cout << "    R² = " << std::fixed << std::setprecision(6) << result.r2 << std::endl;
    std::cout << "    Order = " << result.order << std::endl;
    std::cout << "    Model: " << model.num_exp_modes() << " exp + " 
              << model.num_osc_modes() << " osc modes" << std::endl;
    std::cout << "    Max absolute error = " << max_error << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Zero kernel (edge case)
// =============================================================================
bool TestZeroKernel() {
    std::cout << "  Testing zero kernel (edge case)..." << std::endl;

    const int N = 100;
    Eigen::VectorXd K = Eigen::VectorXd::Zero(N);

    StateSpaceOptions opts;
    opts.max_order = 5;
    opts.r2_threshold = 0.95;

    RadiationStateSpaceFitter fitter(opts);
    StateSpaceFitResult result = fitter.FitKernel(K, 0.01);

    TEST_ASSERT(result.order == 0, "Zero kernel should return order = 0");
    TEST_ASSERT(!result.IsValid(), "Zero kernel fit should not be valid");

    std::cout << "    Order = " << result.order << " (expected 0)" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Short kernel (edge case)
// =============================================================================
bool TestShortKernel() {
    std::cout << "  Testing short kernel (edge case)..." << std::endl;

    Eigen::VectorXd K(3);
    K << 1.0, 0.5, 0.25;

    StateSpaceOptions opts;
    opts.max_order = 5;
    opts.r2_threshold = 0.95;

    RadiationStateSpaceFitter fitter(opts);
    StateSpaceFitResult result = fitter.FitKernel(K, 0.01);

    TEST_ASSERT(result.order == 0, "Too-short kernel should return order = 0");

    std::cout << "    Order = " << result.order << " (expected 0 for N=3)" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Integration with RadiationStateSpaceModel
// =============================================================================
bool TestIntegrationWithModel() {
    std::cout << "  Testing integration with RadiationStateSpaceModel..." << std::endl;

    const double H_true = 10.0;
    const double alpha_true = 2.0;
    const double dt = 0.01;
    const int N = 500;

    Eigen::VectorXd K = GenerateSingleExpKernel(H_true, alpha_true, dt, N);

    StateSpaceOptions opts;
    opts.max_order = 5;
    opts.r2_threshold = 0.99;

    RadiationStateSpaceFitter fitter(opts);
    StateSpaceFitResult result = fitter.FitKernel(K, dt);

    TEST_ASSERT(result.IsValid(), "Fit should be valid");

    // Create model from fit result
    RadiationStateSpaceModel model = RadiationStateSpaceModel::FromFitResult(result);
    model.Reset();

    // Apply constant velocity and check steady-state force
    const double v0 = 1.0;

    for (int step = 0; step < 1000; ++step) {
        model.Step(dt, v0);
    }

    double force = model.GetForce();

    // For exponential kernel K(t) = H*exp(-α*t), steady-state force with v=v0 is:
    // f_∞ = ∫[0,∞] K(τ)*v0 dτ = v0*H/α
    double expected_force = H_true * v0 / alpha_true;

    TEST_ASSERT_NEAR(force, expected_force, 0.1,
                     "Model steady-state force should match kernel integral");

    std::cout << "    Fitted order = " << result.order << std::endl;
    std::cout << "    Steady-state force = " << force 
              << " (expected: " << expected_force << ")" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: ReconstructKernel and ComputeR2 helpers
// =============================================================================
bool TestHelperFunctions() {
    std::cout << "  Testing helper functions..." << std::endl;

    const double alpha = 2.0;
    const double H_val = 5.0;
    const double dt = 0.1;
    const int N = 50;
    
    Eigen::VectorXd K_original(N);
    for (int k = 0; k < N; ++k) {
        K_original(k) = H_val * std::exp(-alpha * k * dt);
    }
    
    StateSpaceOptions opts;
    opts.max_order = 5;
    opts.r2_threshold = 0.999;
    RadiationStateSpaceFitter fitter(opts);
    StateSpaceFitResult result = fitter.FitKernel(K_original, dt);
    
    Eigen::VectorXd K_reconstructed = 
        RadiationStateSpaceFitter::ReconstructKernel(result, dt, N);

    TEST_ASSERT_NEAR(K_reconstructed(0), H_val, 0.1, "K[0] should be ~5.0");
    
    double expected_K1 = H_val * std::exp(-alpha * dt);
    TEST_ASSERT_NEAR(K_reconstructed(1), expected_K1, 0.1, "K[1] should match");

    // Test ComputeR2
    Eigen::VectorXd K_actual(4);
    K_actual << 1.0, 2.0, 3.0, 4.0;
    
    double r2_perfect = RadiationStateSpaceFitter::ComputeR2(K_actual, K_actual);
    TEST_ASSERT_NEAR(r2_perfect, 1.0, 1e-10, "R² of perfect fit should be 1.0");

    double mean = K_actual.mean();
    Eigen::VectorXd K_mean = Eigen::VectorXd::Constant(4, mean);
    double r2_mean = RadiationStateSpaceFitter::ComputeR2(K_actual, K_mean);
    TEST_ASSERT_NEAR(r2_mean, 0.0, 1e-10, "R² of mean predictor should be 0.0");

    std::cout << "    ReconstructKernel: K[0]=" << K_reconstructed(0) << ", K[1]=" << K_reconstructed(1) << std::endl;
    std::cout << "    ComputeR2 (perfect): " << r2_perfect << std::endl;
    std::cout << "    ComputeR2 (mean): " << r2_mean << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: R2 threshold stopping
// =============================================================================
bool TestR2ThresholdStopping() {
    std::cout << "  Testing R2 threshold stopping behavior..." << std::endl;

    const double H_true = 5.0;
    const double alpha_true = 2.0;
    const double dt = 0.01;
    const int N = 500;

    Eigen::VectorXd K = GenerateSingleExpKernel(H_true, alpha_true, dt, N);

    // Test with low R² threshold
    StateSpaceOptions opts_low;
    opts_low.max_order = 10;
    opts_low.r2_threshold = 0.90;

    RadiationStateSpaceFitter fitter_low(opts_low);
    StateSpaceFitResult result_low = fitter_low.FitKernel(K, dt);

    // Test with high R² threshold
    StateSpaceOptions opts_high;
    opts_high.max_order = 10;
    opts_high.r2_threshold = 0.999;

    RadiationStateSpaceFitter fitter_high(opts_high);
    StateSpaceFitResult result_high = fitter_high.FitKernel(K, dt);

    TEST_ASSERT(result_low.IsValid(), "Low threshold fit should be valid");
    TEST_ASSERT(result_high.IsValid(), "High threshold fit should be valid");
    TEST_ASSERT(result_low.r2 >= 0.90, "Low threshold R² should be >= 0.90");
    TEST_ASSERT(result_high.r2 >= 0.99, "High threshold R² should be >= 0.99");

    std::cout << "    Low threshold (0.90): order=" << result_low.order 
              << ", R²=" << result_low.r2 << std::endl;
    std::cout << "    High threshold (0.999): order=" << result_high.order 
              << ", R²=" << result_high.r2 << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Main
// =============================================================================
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "RadiationStateSpaceFitter Unit Tests" << std::endl;
    std::cout << "========================================" << std::endl;

    int passed = 0;
    int failed = 0;

    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "\n[TEST] " << name << std::endl;
        if (test_func()) {
            ++passed;
        } else {
            ++failed;
        }
    };

    run_test("SingleExponential", TestSingleExponential);
    run_test("TwoModeExponential", TestTwoModeExponential);
    run_test("OscillatoryKernel", TestOscillatoryKernel);
    run_test("ZeroKernel", TestZeroKernel);
    run_test("ShortKernel", TestShortKernel);
    run_test("IntegrationWithModel", TestIntegrationWithModel);
    run_test("HelperFunctions", TestHelperFunctions);
    run_test("R2ThresholdStopping", TestR2ThresholdStopping);

    std::cout << "\n========================================" << std::endl;
    std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return (failed > 0) ? 1 : 0;
}
