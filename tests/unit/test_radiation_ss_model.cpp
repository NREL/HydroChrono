/*********************************************************************
 * @file  test_radiation_ss_model.cpp
 * @brief Unit tests for RadiationStateSpaceModel.
 *
 * These tests verify the mathematical correctness of the state-space
 * radiation model using simple analytical cases.
 *
 * TEST SCENARIO (Single DOF, Single Mode):
 *   - num_dofs = 1, num_modes = 1
 *   - Decay rate α > 0, gain H > 0
 *   - Constant velocity v(t) = v0
 *
 * EXPECTED BEHAVIOR:
 *   The ODE: ż = v0 - α z
 *   has analytical solution: z(t) = (v0/α) * (1 - exp(-α t))
 *
 *   Therefore:
 *   - Steady-state (t → ∞): z_∞ = v0/α, f_rad_∞ = H * v0 / α
 *   - Transient: z(t) approaches z_∞ exponentially with time constant 1/α
 *********************************************************************/

#include "../../src/hydro/radiation/radiation_ss_model.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
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
// Test: Constructor validation
// =============================================================================
bool TestConstructorValidation() {
    std::cout << "  Testing constructor validation..." << std::endl;

    // Valid construction
    {
        RadiationStateSpaceModel model(6, 3);
        TEST_ASSERT(model.num_dofs() == 6, "num_dofs should be 6");
        TEST_ASSERT(model.num_modes() == 3, "num_modes should be 3");
    }

    // Invalid: zero DOFs
    {
        bool threw = false;
        try {
            RadiationStateSpaceModel model(0, 3);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        TEST_ASSERT(threw, "Should throw for num_dofs = 0");
    }

    // Invalid: negative modes
    {
        bool threw = false;
        try {
            RadiationStateSpaceModel model(6, -1);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        TEST_ASSERT(threw, "Should throw for num_modes < 0");
    }

    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: SetModeParameters validation
// =============================================================================
bool TestSetModeParametersValidation() {
    std::cout << "  Testing SetModeParameters validation..." << std::endl;

    RadiationStateSpaceModel model(2, 2);

    // Valid parameters
    {
        Eigen::VectorXd h(2);
        h << 1.0, 2.0;
        model.SetModeParameters(0, 0.5, h);
        model.SetModeParameters(1, 1.0, h);
        TEST_ASSERT(model.alphas()(0) == 0.5, "alpha[0] should be 0.5");
        TEST_ASSERT(model.alphas()(1) == 1.0, "alpha[1] should be 1.0");
    }

    // Invalid: mode index out of range
    {
        Eigen::VectorXd h(2);
        h << 1.0, 2.0;
        bool threw = false;
        try {
            model.SetModeParameters(2, 0.5, h);  // index 2 is out of range
        } catch (const std::out_of_range&) {
            threw = true;
        }
        TEST_ASSERT(threw, "Should throw for mode_index >= num_modes");
    }

    // Invalid: alpha <= 0
    {
        Eigen::VectorXd h(2);
        h << 1.0, 2.0;
        bool threw = false;
        try {
            model.SetModeParameters(0, 0.0, h);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        TEST_ASSERT(threw, "Should throw for alpha = 0");
    }

    // Invalid: h_column size mismatch
    {
        Eigen::VectorXd h(3);  // Wrong size
        h << 1.0, 2.0, 3.0;
        bool threw = false;
        try {
            model.SetModeParameters(0, 0.5, h);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        TEST_ASSERT(threw, "Should throw for h_column size mismatch");
    }

    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Single DOF, single mode - steady state convergence
// =============================================================================
bool TestSingleDofSingleModeSteadyState() {
    std::cout << "  Testing single DOF, single mode steady state..." << std::endl;

    // Parameters:
    //   α = 2.0   (decay rate, time constant τ = 0.5 s)
    //   H = 10.0  (gain)
    //   v0 = 1.0  (constant velocity)
    //
    // Expected steady-state:
    //   z_∞ = v0 / α = 0.5
    //   f_rad_∞ = H * z_∞ = H * v0 / α = 5.0

    const double alpha = 2.0;
    const double H_val = 10.0;
    const double v0 = 1.0;
    const double expected_z_ss = v0 / alpha;           // 0.5
    const double expected_f_ss = H_val * v0 / alpha;   // 5.0

    RadiationStateSpaceModel model(1, 1);

    Eigen::VectorXd h(1);
    h << H_val;
    model.SetModeParameters(0, alpha, h);
    model.Reset();

    // Simulate with constant velocity for long enough to reach steady state
    // Time constant τ = 1/α = 0.5 s
    // After 5τ = 2.5 s, we should be within 1% of steady state
    const double dt = 0.01;
    const double total_time = 5.0;  // 10 time constants
    const int num_steps = static_cast<int>(total_time / dt);

    Eigen::VectorXd v(1);
    v << v0;

    for (int step = 0; step < num_steps; ++step) {
        model.Step(dt, v);
    }

    // Check steady state
    Eigen::VectorXd forces = model.GetForces();
    
    const double tolerance = 0.01;  // 1% tolerance
    TEST_ASSERT_NEAR(forces(0), expected_f_ss, expected_f_ss * tolerance,
                     "Radiation force should converge to steady state");

    // Also check internal state
    TEST_ASSERT_NEAR(model.z()(0, 0), expected_z_ss, expected_z_ss * tolerance,
                     "Internal state z should converge to v0/alpha");

    std::cout << "    Steady-state force: " << forces(0) 
              << " (expected: " << expected_f_ss << ")" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Single DOF, single mode - transient response
// =============================================================================
bool TestSingleDofSingleModeTransient() {
    std::cout << "  Testing single DOF, single mode transient..." << std::endl;

    // Parameters (same as steady-state test)
    const double alpha = 2.0;
    const double H_val = 10.0;
    const double v0 = 1.0;

    // Analytical solution for z(t) with z(0) = 0 and constant v = v0:
    //   z(t) = (v0/α) * (1 - exp(-α t))
    // Therefore:
    //   f_rad(t) = H * z(t) = (H * v0 / α) * (1 - exp(-α t))

    RadiationStateSpaceModel model(1, 1);

    Eigen::VectorXd h(1);
    h << H_val;
    model.SetModeParameters(0, alpha, h);
    model.Reset();

    const double dt = 0.01;
    Eigen::VectorXd v(1);
    v << v0;

    // Test at specific times
    const double test_times[] = {0.25, 0.5, 1.0, 2.0};
    const int num_test_times = 4;

    double current_time = 0.0;
    int test_index = 0;

    while (test_index < num_test_times) {
        double target_time = test_times[test_index];
        
        // Step until we reach the target time
        while (current_time < target_time - dt/2) {
            model.Step(dt, v);
            current_time += dt;
        }

        // Compute expected value from analytical solution
        double expected_f = (H_val * v0 / alpha) * (1.0 - std::exp(-alpha * current_time));
        double actual_f = model.GetForces()(0);

        // Allow 2% tolerance for numerical integration
        const double tolerance = 0.02;
        TEST_ASSERT_NEAR(actual_f, expected_f, std::abs(expected_f) * tolerance + 1e-10,
                         std::string("Transient at t=") + std::to_string(current_time));

        std::cout << "    t=" << std::fixed << std::setprecision(2) << current_time 
                  << "s: f_rad=" << std::setprecision(4) << actual_f 
                  << " (expected: " << expected_f << ")" << std::endl;

        ++test_index;
    }

    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Multiple modes superposition
// =============================================================================
bool TestMultipleModes() {
    std::cout << "  Testing multiple modes superposition..." << std::endl;

    // Two modes with different decay rates
    // Mode 0: α₀ = 1.0, H₀ = 2.0
    // Mode 1: α₁ = 4.0, H₁ = 3.0
    // Constant velocity v0 = 1.0
    //
    // Steady-state contributions:
    //   f₀_∞ = H₀ * v0 / α₀ = 2.0
    //   f₁_∞ = H₁ * v0 / α₁ = 0.75
    //   f_total_∞ = 2.75

    const double alpha0 = 1.0, H0 = 2.0;
    const double alpha1 = 4.0, H1 = 3.0;
    const double v0 = 1.0;
    const double expected_f_ss = H0 * v0 / alpha0 + H1 * v0 / alpha1;  // 2.75

    RadiationStateSpaceModel model(1, 2);

    Eigen::VectorXd h0(1), h1(1);
    h0 << H0;
    h1 << H1;
    model.SetModeParameters(0, alpha0, h0);
    model.SetModeParameters(1, alpha1, h1);
    model.Reset();

    // Run to steady state
    const double dt = 0.01;
    const int num_steps = 1000;  // 10 seconds, plenty for both modes to settle

    Eigen::VectorXd v(1);
    v << v0;

    for (int step = 0; step < num_steps; ++step) {
        model.Step(dt, v);
    }

    double actual_f = model.GetForces()(0);
    const double tolerance = 0.01;

    TEST_ASSERT_NEAR(actual_f, expected_f_ss, expected_f_ss * tolerance,
                     "Total force should be sum of mode contributions");

    std::cout << "    Total steady-state force: " << actual_f 
              << " (expected: " << expected_f_ss << ")" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Multi-DOF system
// =============================================================================
bool TestMultiDof() {
    std::cout << "  Testing multi-DOF system..." << std::endl;

    // 2 DOFs, 1 mode
    // α = 2.0
    // H = [10.0, 20.0]  (different gains for each DOF)
    // v = [1.0, 0.5]    (different velocities)
    //
    // Expected steady-state:
    //   z₀_∞ = v₀ / α = 0.5
    //   z₁_∞ = v₁ / α = 0.25
    //   f₀_∞ = H₀ * z₀_∞ = 10.0 * 0.5 = 5.0
    //   f₁_∞ = H₁ * z₁_∞ = 20.0 * 0.25 = 5.0

    const double alpha = 2.0;

    RadiationStateSpaceModel model(2, 1);

    Eigen::VectorXd h(2);
    h << 10.0, 20.0;
    model.SetModeParameters(0, alpha, h);
    model.Reset();

    // Run to steady state
    const double dt = 0.01;
    const int num_steps = 500;

    Eigen::VectorXd v(2);
    v << 1.0, 0.5;

    for (int step = 0; step < num_steps; ++step) {
        model.Step(dt, v);
    }

    Eigen::VectorXd forces = model.GetForces();
    const double tolerance = 0.01;

    TEST_ASSERT_NEAR(forces(0), 5.0, 5.0 * tolerance, "DOF 0 force");
    TEST_ASSERT_NEAR(forces(1), 5.0, 5.0 * tolerance, "DOF 1 force");

    std::cout << "    Forces: [" << forces(0) << ", " << forces(1) << "]"
              << " (expected: [5.0, 5.0])" << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

// =============================================================================
// Test: Reset functionality
// =============================================================================
bool TestReset() {
    std::cout << "  Testing reset functionality..." << std::endl;

    RadiationStateSpaceModel model(1, 1);

    Eigen::VectorXd h(1);
    h << 10.0;
    model.SetModeParameters(0, 2.0, h);

    // Run some steps
    Eigen::VectorXd v(1);
    v << 1.0;
    for (int i = 0; i < 100; ++i) {
        model.Step(0.01, v);
    }

    // Force should be non-zero
    double force_before_reset = model.GetForces()(0);
    TEST_ASSERT(std::abs(force_before_reset) > 0.1, "Force should be non-zero before reset");

    // Reset
    model.Reset();

    // Force should be zero
    double force_after_reset = model.GetForces()(0);
    TEST_ASSERT_NEAR(force_after_reset, 0.0, 1e-10, "Force should be zero after reset");

    std::cout << "    Force before reset: " << force_before_reset << std::endl;
    std::cout << "    Force after reset: " << force_after_reset << std::endl;
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

    run_test("ConstructorValidation", TestConstructorValidation);
    run_test("SetModeParametersValidation", TestSetModeParametersValidation);
    run_test("SingleDofSingleModeSteadyState", TestSingleDofSingleModeSteadyState);
    run_test("SingleDofSingleModeTransient", TestSingleDofSingleModeTransient);
    run_test("MultipleModes", TestMultipleModes);
    run_test("MultiDof", TestMultiDof);
    run_test("Reset", TestReset);

    std::cout << "\n========================================" << std::endl;
    std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return (failed > 0) ? 1 : 0;
}

