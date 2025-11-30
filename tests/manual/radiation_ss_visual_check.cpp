/*********************************************************************
 * @file  radiation_ss_visual_check.cpp
 * @brief Diagnostic tool for visual verification of RadiationStateSpaceModel.
 *
 * This program outputs CSV data comparing the numerical model against
 * analytical solutions for simple test cases. Use with the accompanying
 * Python script (plot_radiation_ss.py) for visualization.
 *
 * NOT a unit test - this is for manual/visual verification only.
 *
 * USAGE:
 *   radiation_ss_visual_check [output_dir]
 *   
 *   If output_dir is not specified, writes to ./output/ subdirectory.
 *
 * OUTPUT FILES (in output/ subdirectory):
 *   - radiation_ss_constant_v.csv   : Constant velocity response
 *   - radiation_ss_step_v.csv       : Step velocity response
 *   - radiation_ss_multi_mode.csv   : Multiple modes superposition
 *   - radiation_ss_dt_sensitivity.csv : Time step sensitivity
 *********************************************************************/

#include "../../src/hydro/radiation/radiation_ss_model.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <filesystem>

using namespace hydrochrono::hydro;

/**
 * @brief Scenario 1: Constant velocity response.
 *
 * Parameters:
 *   α = 0.5, H = 2.0, v₀ = 1.0
 *
 * Analytical solution:
 *   z(t) = (v₀/α) * (1 - exp(-α t))
 *   f(t) = H * z(t)
 *
 * Steady state (t → ∞):
 *   z_∞ = v₀/α = 2.0
 *   f_∞ = H * v₀/α = 4.0
 */
void WriteConstantVelocityResponse(const std::string& output_dir) {
    const std::string filename = output_dir + "/radiation_ss_constant_v.csv";
    std::cout << "Writing: " << filename << std::endl;

    const int num_dofs = 1;
    const int num_modes = 1;
    const double alpha = 0.5;   // Time constant τ = 2.0 s
    const double H = 2.0;
    const double v0 = 1.0;
    const double dt = 0.02;
    const int steps = 500;      // 10 seconds total

    RadiationStateSpaceModel model(num_dofs, num_modes);
    model.SetModeParameters(0, alpha, Eigen::VectorXd::Constant(1, H));
    model.Reset();

    std::ofstream file(filename);
    file << "t,v,z_model,z_exact,f_model,f_exact,error_pct\n";
    file << std::fixed << std::setprecision(6);

    Eigen::VectorXd v(1);
    v(0) = v0;

    double t = 0.0;

    // Write initial state (before first step)
    file << t << "," << v0 << ",0,0,0,0,0\n";

    for (int k = 0; k < steps; ++k) {
        model.Step(dt, v);
        t += dt;

        double z_model = model.z()(0, 0);
        double f_model = model.GetForces()(0);

        double z_exact = (v0 / alpha) * (1.0 - std::exp(-alpha * t));
        double f_exact = H * z_exact;

        double error_pct = (std::abs(f_exact) > 1e-10) 
            ? 100.0 * std::abs(f_model - f_exact) / std::abs(f_exact) 
            : 0.0;

        file << t << "," << v0 << "," << z_model << "," << z_exact << ","
             << f_model << "," << f_exact << "," << error_pct << "\n";
    }

    std::cout << "  -> Steady state: f_∞ = " << H * v0 / alpha << std::endl;
}

/**
 * @brief Scenario 2: Step velocity response (velocity changes at t=5s).
 *
 * v(t) = v₀ for t < 5s, then v(t) = v₁ for t >= 5s
 *
 * This tests the model's response to a sudden velocity change.
 */
void WriteStepVelocityResponse(const std::string& output_dir) {
    const std::string filename = output_dir + "/radiation_ss_step_v.csv";
    std::cout << "Writing: " << filename << std::endl;

    const int num_dofs = 1;
    const int num_modes = 1;
    const double alpha = 1.0;   // Time constant τ = 1.0 s
    const double H = 5.0;
    const double v0 = 1.0;
    const double v1 = 0.0;      // Step to zero at t=5s
    const double t_step = 5.0;
    const double dt = 0.02;
    const int steps = 500;      // 10 seconds total

    RadiationStateSpaceModel model(num_dofs, num_modes);
    model.SetModeParameters(0, alpha, Eigen::VectorXd::Constant(1, H));
    model.Reset();

    std::ofstream file(filename);
    file << "t,v,z_model,f_model\n";
    file << std::fixed << std::setprecision(6);

    Eigen::VectorXd v(1);
    double t = 0.0;

    file << t << "," << v0 << ",0,0\n";

    for (int k = 0; k < steps; ++k) {
        // Switch velocity at t_step
        v(0) = (t < t_step) ? v0 : v1;

        model.Step(dt, v);
        t += dt;

        double z_model = model.z()(0, 0);
        double f_model = model.GetForces()(0);

        file << t << "," << v(0) << "," << z_model << "," << f_model << "\n";
    }

    std::cout << "  -> Step from v=" << v0 << " to v=" << v1 << " at t=" << t_step << "s" << std::endl;
}

/**
 * @brief Scenario 3: Multiple modes superposition.
 *
 * Two modes with different decay rates:
 *   Mode 0: α₀ = 0.5, H₀ = 3.0  (slow mode, τ = 2s)
 *   Mode 1: α₁ = 2.0, H₁ = 1.0  (fast mode, τ = 0.5s)
 *
 * The fast mode contributes quickly but saturates at a lower level.
 * The slow mode takes longer but dominates at steady state.
 */
void WriteMultipleModeResponse(const std::string& output_dir) {
    const std::string filename = output_dir + "/radiation_ss_multi_mode.csv";
    std::cout << "Writing: " << filename << std::endl;

    const int num_dofs = 1;
    const int num_modes = 2;
    
    const double alpha0 = 0.5, H0 = 3.0;  // Slow mode
    const double alpha1 = 2.0, H1 = 1.0;  // Fast mode
    const double v0 = 1.0;
    const double dt = 0.02;
    const int steps = 600;  // 12 seconds

    RadiationStateSpaceModel model(num_dofs, num_modes);
    
    Eigen::VectorXd h0(1), h1(1);
    h0 << H0;
    h1 << H1;
    model.SetModeParameters(0, alpha0, h0);
    model.SetModeParameters(1, alpha1, h1);
    model.Reset();

    std::ofstream file(filename);
    file << "t,v,z0_model,z1_model,f_mode0,f_mode1,f_total,f0_exact,f1_exact,f_total_exact\n";
    file << std::fixed << std::setprecision(6);

    Eigen::VectorXd v(1);
    v(0) = v0;

    double t = 0.0;
    file << t << "," << v0 << ",0,0,0,0,0,0,0,0\n";

    for (int k = 0; k < steps; ++k) {
        model.Step(dt, v);
        t += dt;

        double z0 = model.z()(0, 0);
        double z1 = model.z()(0, 1);
        double f_total = model.GetForces()(0);

        // Individual mode contributions
        double f_mode0 = H0 * z0;
        double f_mode1 = H1 * z1;

        // Analytical solutions
        double z0_exact = (v0 / alpha0) * (1.0 - std::exp(-alpha0 * t));
        double z1_exact = (v0 / alpha1) * (1.0 - std::exp(-alpha1 * t));
        double f0_exact = H0 * z0_exact;
        double f1_exact = H1 * z1_exact;
        double f_total_exact = f0_exact + f1_exact;

        file << t << "," << v0 << "," << z0 << "," << z1 << ","
             << f_mode0 << "," << f_mode1 << "," << f_total << ","
             << f0_exact << "," << f1_exact << "," << f_total_exact << "\n";
    }

    double f_ss_0 = H0 * v0 / alpha0;
    double f_ss_1 = H1 * v0 / alpha1;
    std::cout << "  -> Mode 0 (slow): τ=" << 1.0/alpha0 << "s, f_∞=" << f_ss_0 << std::endl;
    std::cout << "  -> Mode 1 (fast): τ=" << 1.0/alpha1 << "s, f_∞=" << f_ss_1 << std::endl;
    std::cout << "  -> Total steady-state: f_∞=" << f_ss_0 + f_ss_1 << std::endl;
}

/**
 * @brief Scenario 4: Time step sensitivity analysis.
 *
 * Runs the same scenario with different time steps to verify
 * that the exact exponential integrator remains accurate.
 */
void WriteTimeStepSensitivity(const std::string& output_dir) {
    const std::string filename = output_dir + "/radiation_ss_dt_sensitivity.csv";
    std::cout << "Writing: " << filename << std::endl;

    const int num_dofs = 1;
    const int num_modes = 1;
    const double alpha = 2.0;
    const double H = 5.0;
    const double v0 = 1.0;
    const double total_time = 5.0;

    std::ofstream file(filename);
    file << "t,f_dt001,f_dt01,f_dt05,f_dt1,f_exact\n";
    file << std::fixed << std::setprecision(6);

    // Set up models with different time steps
    double dts[] = {0.01, 0.1, 0.5, 1.0};
    const int num_dts = 4;
    
    RadiationStateSpaceModel models[4] = {
        RadiationStateSpaceModel(num_dofs, num_modes),
        RadiationStateSpaceModel(num_dofs, num_modes),
        RadiationStateSpaceModel(num_dofs, num_modes),
        RadiationStateSpaceModel(num_dofs, num_modes)
    };

    Eigen::VectorXd h(1);
    h << H;
    for (int i = 0; i < num_dts; ++i) {
        models[i].SetModeParameters(0, alpha, h);
        models[i].Reset();
    }

    Eigen::VectorXd v(1);
    v(0) = v0;

    // Sample at 0.01s intervals
    const double dt_sample = 0.01;
    double t = 0.0;
    double model_times[4] = {0, 0, 0, 0};

    file << "0,0,0,0,0,0\n";

    int total_samples = static_cast<int>(total_time / dt_sample);
    for (int sample = 0; sample < total_samples; ++sample) {
        t += dt_sample;

        // Advance each model to time t
        double forces[4];
        for (int i = 0; i < num_dts; ++i) {
            while (model_times[i] + dts[i] <= t + 1e-9) {
                models[i].Step(dts[i], v);
                model_times[i] += dts[i];
            }
            forces[i] = models[i].GetForces()(0);
        }

        double f_exact = (H * v0 / alpha) * (1.0 - std::exp(-alpha * t));

        file << t << "," << forces[0] << "," << forces[1] << ","
             << forces[2] << "," << forces[3] << "," << f_exact << "\n";
    }

    std::cout << "  -> Comparing dt = 0.01, 0.1, 0.5, 1.0 s" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string output_dir = "output";
    if (argc > 1) {
        output_dir = argv[1];
    }

    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_dir);

    std::cout << "======================================================" << std::endl;
    std::cout << "RadiationStateSpaceModel Visual Verification" << std::endl;
    std::cout << "======================================================" << std::endl;
    std::cout << "Output directory: " << std::filesystem::absolute(output_dir) << std::endl;
    std::cout << std::endl;

    WriteConstantVelocityResponse(output_dir);
    std::cout << std::endl;

    WriteStepVelocityResponse(output_dir);
    std::cout << std::endl;

    WriteMultipleModeResponse(output_dir);
    std::cout << std::endl;

    WriteTimeStepSensitivity(output_dir);
    std::cout << std::endl;

    std::cout << "======================================================" << std::endl;
    std::cout << "Done! Use plot_radiation_ss.py to visualize results." << std::endl;
    std::cout << "======================================================" << std::endl;

    return 0;
}

