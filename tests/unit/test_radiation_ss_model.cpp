/*********************************************************************
 * @file  test_radiation_ss_model.cpp
 * @brief Unit tests for RadiationStateSpaceModel.
 *
 * These tests verify the mathematical correctness of the state-space
 * radiation model for both exponential and oscillatory modes.
 *
 * TEST SCENARIOS:
 *   1. Exponential mode with constant velocity
 *   2. Oscillatory mode with constant velocity  
 *   3. Multiple modes combined
 *   4. FromFitResult factory method
 *   5. ReconstructKernel verification
 *
 * ANALYTICAL SOLUTIONS:
 *   Exponential mode (z' = -α*z + b*v):
 *     z(t) = (b*v/α) * (1 - exp(-α*t))
 *     z_∞ = b*v/α
 *     f_∞ = H * b * v / α
 *
 *   Oscillatory mode (damped oscillation):
 *     More complex, verified via kernel reconstruction
 *********************************************************************/

#include "../../src/hydro/radiation/radiation_ss_model.h"
#include "../../src/hydro/radiation/radiation_ss_fitter.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

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
// Test: Exponential mode steady-state
// =============================================================================
bool TestExponentialModeSteadyState() {
    std::cout << "  Testing exponential mode steady-state..." << std::endl;

    // Parameters
    const double alpha = 2.0;  // Decay rate
    const double b = 1.0;      // Input gain
    const double H = 10.0;     // Output gain
    const double v = 1.0;      // Constant velocity
    const double dt = 0.01;
    const int num_steps = 1000;

    // Expected steady-state: z_∞ = b*v/α, f_∞ = H*b*v/α
    const double z_inf = b * v / alpha;
    const double f_inf = H * z_inf;

    // Create model
    RadiationStateSpaceModel model;
    model.AddExponentialMode(alpha, b, H);
    model.Reset();

    // Run to steady state
    for (int step = 0; step < num_steps; ++step) {
        model.Step(dt, v);
    }

    double force = model.GetForce();
    const double tolerance = 0.01 * std::abs(f_inf);

    TEST_ASSERT_NEAR(force, f_inf, tolerance, "Steady-state force");

    std::cout << "    Force: " << force << " (expected: " << f_inf << ")" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Exponential mode transient
// =============================================================================
bool TestExponentialModeTransient() {
    std::cout << "  Testing exponential mode transient..." << std::endl;

    const double alpha = 2.0;
    const double b = 1.0;
    const double H = 10.0;
    const double v = 1.0;
    const double dt = 0.001;  // Small dt for accuracy

    RadiationStateSpaceModel model;
    model.AddExponentialMode(alpha, b, H);
    model.Reset();

    // Test at specific times
    std::vector<double> test_times = {0.5, 1.0, 2.0, 3.0};
    
    for (double t_target : test_times) {
        model.Reset();
        int num_steps = static_cast<int>(t_target / dt);
        
        for (int step = 0; step < num_steps; ++step) {
            model.Step(dt, v);
        }

        // Expected: z(t) = (b*v/α) * (1 - exp(-α*t))
        // f(t) = H * z(t) = H*b*v/α * (1 - exp(-α*t))
        double z_exact = (b * v / alpha) * (1.0 - std::exp(-alpha * t_target));
        double f_exact = H * z_exact;
        double f_model = model.GetForce();

        // Allow 1% tolerance
        double tolerance = 0.01 * std::abs(f_exact) + 1e-10;
        TEST_ASSERT_NEAR(f_model, f_exact, tolerance, 
            "Force at t=" + std::to_string(t_target));

        std::cout << "    t=" << t_target << ": f=" << f_model 
                  << " (exact: " << f_exact << ")" << std::endl;
    }

    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Oscillatory mode kernel reconstruction
// =============================================================================
bool TestOscillatoryModeKernel() {
    std::cout << "  Testing oscillatory mode kernel reconstruction..." << std::endl;

    const double alpha = 0.8;   // Decay rate
    const double omega = 3.0;   // Oscillation frequency
    const double b_c = 1.0;     // Input gain (cos)
    const double b_s = 0.0;     // Input gain (sin)
    const double H_c = 5.0;     // Output gain (cos)
    const double H_s = 0.0;     // Output gain (sin)
    const double dt = 0.02;
    const int num_samples = 100;

    RadiationStateSpaceModel model;
    model.AddOscillatoryMode(alpha, omega, b_c, b_s, H_c, H_s);

    // Reconstruct kernel
    Eigen::VectorXd K = model.ReconstructKernel(dt, num_samples);

    // Expected kernel: K(t) = H_c * exp(-α*t) * cos(ω*t)
    // (since b_c=1, b_s=0, H_s=0)
    double max_error = 0.0;
    for (int k = 0; k < num_samples; ++k) {
        double t = k * dt;
        double K_exact = H_c * std::exp(-alpha * t) * std::cos(omega * t);
        double error = std::abs(K(k) - K_exact);
        max_error = std::max(max_error, error);
    }

    const double tolerance = 1e-10;
    TEST_ASSERT(max_error < tolerance, 
        "Kernel reconstruction error: " + std::to_string(max_error));

    std::cout << "    Max kernel error: " << max_error << std::endl;
    std::cout << "    K(0) = " << K(0) << " (expected: " << H_c << ")" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Oscillatory mode force response
// =============================================================================
bool TestOscillatoryModeForce() {
    std::cout << "  Testing oscillatory mode force response..." << std::endl;

    const double alpha = 0.8;
    const double omega = 3.0;
    const double H_c = 5.0;
    const double H_s = 0.0;
    const double dt = 0.01;
    const int num_steps = 500;

    RadiationStateSpaceModel model;
    model.AddOscillatoryMode(alpha, omega, 1.0, 0.0, H_c, H_s);
    model.Reset();

    // Apply impulse (velocity = 1 for one step, then 0)
    // This should show clear oscillatory decay in the force
    model.Step(dt, 1.0);
    
    double prev_force = model.GetForce();
    int sign_changes = 0;
    double max_force = std::abs(prev_force);

    for (int step = 1; step < num_steps; ++step) {
        model.Step(dt, 0.0);  // No input after initial impulse
        double force = model.GetForce();
        
        if (prev_force * force < 0) {
            ++sign_changes;
        }
        prev_force = force;
        max_force = std::max(max_force, std::abs(force));
    }

    // With omega=3, over 5 seconds we expect ~2.4 full oscillations = ~4-5 sign changes
    TEST_ASSERT(sign_changes >= 3, "Force should oscillate (at least 3 sign changes)");
    TEST_ASSERT(max_force > 0.01, "Force should have significant amplitude");

    std::cout << "    Sign changes: " << sign_changes << std::endl;
    std::cout << "    Max force: " << max_force << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Multiple modes
// =============================================================================
bool TestMultipleModes() {
    std::cout << "  Testing multiple modes..." << std::endl;

    RadiationStateSpaceModel model;
    
    // Add two exponential modes
    model.AddExponentialMode(1.0, 1.0, 2.0);  // α=1, b=1, H=2
    model.AddExponentialMode(2.0, 1.0, 3.0);  // α=2, b=1, H=3
    model.Reset();

    const double v = 1.0;
    const double dt = 0.01;
    const int num_steps = 1000;

    // Run to steady state
    for (int step = 0; step < num_steps; ++step) {
        model.Step(dt, v);
    }

    // Expected steady-state force:
    // f = H1*b1*v/α1 + H2*b2*v/α2 = 2*1*1/1 + 3*1*1/2 = 2 + 1.5 = 3.5
    double f_expected = 2.0 + 1.5;
    double f_model = model.GetForce();

    const double tolerance = 0.01 * f_expected;
    TEST_ASSERT_NEAR(f_model, f_expected, tolerance, "Multiple mode steady-state");

    std::cout << "    Force: " << f_model << " (expected: " << f_expected << ")" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Mixed exponential and oscillatory modes
// =============================================================================
bool TestMixedModes() {
    std::cout << "  Testing mixed exponential + oscillatory modes..." << std::endl;

    RadiationStateSpaceModel model;
    model.AddExponentialMode(1.0, 1.0, 5.0);        // Exp mode
    model.AddOscillatoryMode(0.5, 2.0, 1.0, 0.0, 3.0, 0.0);  // Osc mode
    model.Reset();

    // Verify kernel at t=0
    Eigen::VectorXd K = model.ReconstructKernel(0.01, 10);
    
    // At t=0: K(0) = H_exp*b_exp + H_c_osc*b_c_osc = 5*1 + 3*1 = 8
    double K_0_expected = 5.0 + 3.0;
    TEST_ASSERT_NEAR(K(0), K_0_expected, 1e-10, "Mixed mode K(0)");

    std::cout << "    K(0) = " << K(0) << " (expected: " << K_0_expected << ")" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Reset functionality
// =============================================================================
bool TestReset() {
    std::cout << "  Testing reset functionality..." << std::endl;

    RadiationStateSpaceModel model;
    model.AddExponentialMode(2.0, 1.0, 10.0);

    // Run some steps
    for (int i = 0; i < 100; ++i) {
        model.Step(0.01, 1.0);
    }

    // Force should be non-zero
    double force_before = model.GetForce();
    TEST_ASSERT(std::abs(force_before) > 0.1, "Force should be non-zero before reset");

    // Reset
    model.Reset();

    // Force should be zero
    double force_after = model.GetForce();
    TEST_ASSERT_NEAR(force_after, 0.0, 1e-10, "Force should be zero after reset");

    std::cout << "    Force before reset: " << force_before << std::endl;
    std::cout << "    Force after reset: " << force_after << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: FromFitResult factory
// =============================================================================
bool TestFromFitResult() {
    std::cout << "  Testing FromFitResult factory..." << std::endl;

    // Create a simple 2x2 state-space system with oscillatory eigenvalues
    StateSpaceFitResult result;
    result.order = 2;
    result.A.resize(2, 2);
    result.B.resize(2);
    result.C.resize(1, 2);

    // A matrix with eigenvalues -0.5 ± 2j (oscillatory decay)
    double alpha = 0.5;
    double omega = 2.0;
    result.A << -alpha, -omega,
                 omega, -alpha;
    result.B << 1.0, 0.0;
    result.C << 1.0, 0.0;
    result.D = 0.0;
    result.r2 = 1.0;

    // Create model from fit result
    RadiationStateSpaceModel model = RadiationStateSpaceModel::FromFitResult(result);

    TEST_ASSERT(model.num_osc_modes() == 1, "Should have 1 oscillatory mode");
    TEST_ASSERT(model.num_exp_modes() == 0, "Should have 0 exponential modes");
    TEST_ASSERT(model.total_states() == 2, "Should have 2 total states");

    // Check that kernel reconstruction works
    Eigen::VectorXd K = model.ReconstructKernel(0.01, 100);
    TEST_ASSERT(K(0) > 0, "K(0) should be positive");
    
    std::cout << "    Osc modes: " << model.num_osc_modes() << std::endl;
    std::cout << "    Exp modes: " << model.num_exp_modes() << std::endl;
    std::cout << "    K(0) = " << K(0) << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Invalid parameters throw
// =============================================================================
bool TestInvalidParameters() {
    std::cout << "  Testing invalid parameter validation..." << std::endl;

    RadiationStateSpaceModel model;

    // Test: alpha <= 0 should throw
    {
        bool threw = false;
        try {
            model.AddExponentialMode(0.0, 1.0, 1.0);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        TEST_ASSERT(threw, "Should throw for alpha = 0");
    }

    // Test: omega <= 0 should throw
    {
        bool threw = false;
        try {
            model.AddOscillatoryMode(1.0, 0.0, 1.0, 0.0, 1.0, 0.0);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        TEST_ASSERT(threw, "Should throw for omega = 0");
    }

    // Test: dt <= 0 should throw
    {
        bool threw = false;
        model.AddExponentialMode(1.0, 1.0, 1.0);
        try {
            model.Step(0.0, 1.0);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        TEST_ASSERT(threw, "Should throw for dt = 0");
    }

    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Main
// =============================================================================
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "RadiationStateSpaceModel Unit Tests" << std::endl;
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

    run_test("ExponentialModeSteadyState", TestExponentialModeSteadyState);
    run_test("ExponentialModeTransient", TestExponentialModeTransient);
    run_test("OscillatoryModeKernel", TestOscillatoryModeKernel);
    run_test("OscillatoryModeForce", TestOscillatoryModeForce);
    run_test("MultipleModes", TestMultipleModes);
    run_test("MixedModes", TestMixedModes);
    run_test("Reset", TestReset);
    run_test("FromFitResult", TestFromFitResult);
    run_test("InvalidParameters", TestInvalidParameters);

    std::cout << "\n========================================" << std::endl;
    std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return (failed > 0) ? 1 : 0;
}
