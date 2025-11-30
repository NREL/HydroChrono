/*********************************************************************
 * @file  sphere_decay_ss_comparison.cpp
 * @brief Compare sphere heave decay: Convolution vs State-Space radiation.
 *
 * This diagnostic runs the IEA sphere decay test twice:
 *   1. Using RIRF convolution (baseline)
 *   2. Using state-space approximation
 *
 * Outputs:
 *   - output/sphere_kernel_fit.csv         (t, K_actual, K_fit)
 *   - output/sphere_decay_convolution.csv  (time, heave_z, heave_vz)
 *   - output/sphere_decay_state_space.csv  (time, heave_z, heave_vz)
 *   - output/sphere_performance.csv        (scaling data)
 *
 * Visualize with: python plot_model_comparison.py sphere --save
 *********************************************************************/

#include <hydroc/hydro_system.h>
#include <hydroc/helper.h>
#include <hydroc/io/h5_reader.h>
#include <hydroc/radiation/radiation_types.h>
#include "../src/hydro/radiation/radiation_ss_fitter.h"
#include "../src/hydro/radiation/radiation_ss_model.h"

#include <chrono/core/ChRealtimeStep.h>
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/solver/ChSolver.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace chrono;

// =============================================================================
// Physical Constants for IEA OES Task 10 Sphere
// =============================================================================
namespace {
constexpr double kSphereMass = 261.8e3;           // kg
constexpr double kInitialHeaveDisplacement = -1.0; // m (submerged 1m)
constexpr double kDefaultTimestep = 0.015;        // s
constexpr double kSimulationDuration = 60.0;      // s
}  // namespace

// Ensure directory exists
void ensure_directory_exists(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path);
    }
}

// Extract and fit heave kernel (DOF 2 = heave, self-coupling)
void WriteKernelFit(const std::string& h5_path, const std::filesystem::path& output_dir) {
    std::cout << "\n=== Extracting and Fitting Heave Kernel ===" << std::endl;
    
    // Load H5 file
    H5FileInfo h5_file(h5_path);
    HydroData hydro_data = h5_file.ReadH5Data();
    
    // Get RIRF dimensions and time
    const int rirf_steps = hydro_data.GetRIRFDims(2);
    const Eigen::VectorXd rirf_time = hydro_data.GetRIRFTimeVector();
    const double dt = (rirf_steps > 1) ? (rirf_time(rirf_steps-1) - rirf_time(0)) / (rirf_steps - 1) : 0.01;
    
    // Extract heave-heave kernel (body 0, DOF 2 = heave, column 2)
    const int body = 0;
    const int row_dof = 2;  // heave
    const int col = 2;      // heave
    
    Eigen::VectorXd K(rirf_steps);
    for (int step = 0; step < rirf_steps; ++step) {
        K(step) = hydro_data.GetRIRFVal(body, row_dof, col, step);
    }
    
    std::cout << "  RIRF samples: " << rirf_steps << ", dt = " << dt << " s" << std::endl;
    std::cout << "  K[0] = " << K(0) << ", K[end] = " << K(rirf_steps-1) << std::endl;
    
    // Fit using state-space fitter
    hydrochrono::hydro::StateSpaceOptions opts;
    opts.max_order = 10;
    opts.r2_threshold = 0.99;
    opts.max_hankel_size = 200;
    
    hydrochrono::hydro::RadiationStateSpaceFitter fitter(opts);
    hydrochrono::hydro::StateSpaceFitResult result = fitter.FitKernel(K, dt);
    
    std::cout << "  Fit order: " << result.order << ", R² = " << result.r2 << std::endl;
    
    // Reconstruct kernel
    Eigen::VectorXd K_fit = hydrochrono::hydro::RadiationStateSpaceFitter::ReconstructKernel(result, dt, rirf_steps);
    
    // Write to CSV
    std::ofstream file((output_dir / "sphere_kernel_fit.csv").string());
    file << "t,K_actual,K_fit\n";
    file << std::setprecision(8);
    for (int i = 0; i < rirf_steps; ++i) {
        file << rirf_time(i) << "," << K(i) << "," << K_fit(i) << "\n";
    }
    std::cout << "  Wrote: " << (output_dir / "sphere_kernel_fit.csv").string() << std::endl;
}

// Run performance scaling test
void WritePerformanceScaling(const std::string& h5_path, const std::string& mesh_path,
                             const std::filesystem::path& output_dir) {
    std::cout << "\n=== Running Performance Scaling Test ===" << std::endl;
    
    std::vector<int> step_counts = {100, 500, 1000, 2000, 5000, 10000};
    std::vector<double> conv_times, ss_times;
    
    const double dt = 0.015;
    
    for (int steps : step_counts) {
        
        // Run convolution
        ChSystemNSC system1;
        system1.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
        system1.SetSolverType(ChSolver::Type::GMRES);
        
        auto body1 = chrono_types::make_shared<ChBodyEasyMesh>(mesh_path, 1000, false, true, false);
        body1->SetName("body1");
        body1->SetPos(ChVector3d(0, 0, kInitialHeaveDisplacement));
        body1->SetMass(kSphereMass);
        system1.Add(body1);
        
        std::vector<std::shared_ptr<ChBody>> bodies1 = {body1};
        HydroSystem hydro1(bodies1, h5_path);
        hydro1.SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kRirfConvolution);
        hydro1.AddWaves(std::make_shared<NoWave>(1));
        
        system1.DoStepDynamics(dt * 0.001);  // init
        system1.SetChTime(0.0);
        body1->SetPos(ChVector3d(0, 0, kInitialHeaveDisplacement));
        body1->SetPosDt(ChVector3d(0, 0, 0));
        
        auto t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < steps; ++i) {
            system1.DoStepDynamics(dt);
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        double conv_time = std::chrono::duration<double, std::micro>(t2 - t1).count();
        conv_times.push_back(conv_time);
        
        // Run state-space
        ChSystemNSC system2;
        system2.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
        system2.SetSolverType(ChSolver::Type::GMRES);
        
        auto body2 = chrono_types::make_shared<ChBodyEasyMesh>(mesh_path, 1000, false, true, false);
        body2->SetName("body1");
        body2->SetPos(ChVector3d(0, 0, kInitialHeaveDisplacement));
        body2->SetMass(kSphereMass);
        system2.Add(body2);
        
        std::vector<std::shared_ptr<ChBody>> bodies2 = {body2};
        HydroSystem hydro2(bodies2, h5_path);
        hydro2.SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kStateSpace);
        hydrochrono::hydro::StateSpaceOptions ss_opts;
        ss_opts.max_order = 10;
        ss_opts.r2_threshold = 0.95;
        ss_opts.max_hankel_size = 200;
        hydro2.SetStateSpaceOptions(ss_opts);
        hydro2.AddWaves(std::make_shared<NoWave>(1));
        
        system2.DoStepDynamics(dt * 0.001);  // init (includes fitting)
        system2.SetChTime(0.0);
        body2->SetPos(ChVector3d(0, 0, kInitialHeaveDisplacement));
        body2->SetPosDt(ChVector3d(0, 0, 0));
        
        auto t3 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < steps; ++i) {
            system2.DoStepDynamics(dt);
        }
        auto t4 = std::chrono::high_resolution_clock::now();
        double ss_time = std::chrono::duration<double, std::micro>(t4 - t3).count();
        ss_times.push_back(ss_time);
        
        std::cout << "  Steps=" << steps << ": Conv=" << conv_time/1000 << "ms, SS=" << ss_time/1000 << "ms" << std::endl;
    }
    
    // Write to CSV
    std::ofstream file((output_dir / "sphere_performance.csv").string());
    file << "sim_steps,conv_time_us,ss_time_us,speedup,conv_per_step_us,ss_per_step_us\n";
    file << std::setprecision(6);
    for (size_t i = 0; i < step_counts.size(); ++i) {
        double speedup = conv_times[i] / ss_times[i];
        double conv_per = conv_times[i] / step_counts[i];
        double ss_per = ss_times[i] / step_counts[i];
        file << step_counts[i] << "," << conv_times[i] << "," << ss_times[i] << ","
             << speedup << "," << conv_per << "," << ss_per << "\n";
    }
    std::cout << "  Wrote: " << (output_dir / "sphere_performance.csv").string() << std::endl;
}

// Result storage
struct SimulationResult {
    std::vector<double> time;
    std::vector<double> heave_z;
    std::vector<double> heave_vz;
    double init_time_ms;    // Time for initialization (fitting for SS)
    double sim_time_ms;     // Time for actual simulation loop
    double total_time_ms;   // Total wall time
};

// Run sphere decay simulation with specified radiation method
SimulationResult RunSphereDecay(
    const std::string& h5_path,
    const std::string& mesh_path,
    hydrochrono::hydro::RadiationMethod method,
    const std::string& method_name) {
    
    std::cout << "\n=== Running Sphere Decay: " << method_name << " ===" << std::endl;
    
    SimulationResult result;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // System setup
    ChSystemNSC system;
    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
    
    system.SetSolverType(ChSolver::Type::GMRES);
    system.GetSolver()->AsIterative()->SetMaxIterations(300);
    
    // Create sphere body
    auto sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(
        mesh_path,
        1000,   // density (ignored, mass set manually)
        false,  // don't evaluate mass
        true,   // visualization
        false   // no collision
    );
    
    sphereBody->SetName("body1");
    sphereBody->SetPos(ChVector3d(0, 0, kInitialHeaveDisplacement));
    sphereBody->SetMass(kSphereMass);
    system.Add(sphereBody);
    
    // Setup hydrodynamics - this is where SS fitting happens
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(sphereBody);
    
    auto init_start = std::chrono::high_resolution_clock::now();
    
    HydroSystem hydro_forces(bodies, h5_path);
    
    // Set radiation method
    hydro_forces.SetRadiationMethod(method);
    if (method == hydrochrono::hydro::RadiationMethod::kStateSpace) {
        hydrochrono::hydro::StateSpaceOptions opts;
        opts.max_order = 10;
        opts.r2_threshold = 0.99;
        hydro_forces.SetStateSpaceOptions(opts);
    }
    
    // No waves for decay test
    auto no_waves = std::make_shared<NoWave>(1);
    hydro_forces.AddWaves(no_waves);
    
    // Force initialization of radiation component (lazy init happens here)
    // Do one dummy step to trigger initialization
    system.DoStepDynamics(kDefaultTimestep * 0.001);  // tiny step
    system.SetChTime(0.0);  // reset time
    sphereBody->SetPos(ChVector3d(0, 0, kInitialHeaveDisplacement));  // reset position
    sphereBody->SetPosDt(ChVector3d(0, 0, 0));  // reset velocity
    
    auto init_end = std::chrono::high_resolution_clock::now();
    result.init_time_ms = std::chrono::duration<double, std::milli>(init_end - init_start).count();
    
    // Pre-allocate result vectors
    int expected_steps = static_cast<int>(kSimulationDuration / kDefaultTimestep) + 10;
    result.time.reserve(expected_steps);
    result.heave_z.reserve(expected_steps);
    result.heave_vz.reserve(expected_steps);
    
    // Timing simulation only
    auto sim_start = std::chrono::high_resolution_clock::now();
    
    // Main simulation loop
    while (system.GetChTime() <= kSimulationDuration) {
        system.DoStepDynamics(kDefaultTimestep);
        
        result.time.push_back(system.GetChTime());
        result.heave_z.push_back(sphereBody->GetPos().z());
        result.heave_vz.push_back(sphereBody->GetPosDt().z());
    }
    
    // Timing end
    auto sim_end = std::chrono::high_resolution_clock::now();
    result.sim_time_ms = std::chrono::duration<double, std::milli>(sim_end - sim_start).count();
    result.total_time_ms = std::chrono::duration<double, std::milli>(sim_end - total_start).count();
    
    std::cout << "  Simulation time: " << kSimulationDuration << " s" << std::endl;
    std::cout << "  Init time: " << result.init_time_ms << " ms" << std::endl;
    std::cout << "  Sim loop time: " << result.sim_time_ms << " ms" << std::endl;
    std::cout << "  Total time: " << result.total_time_ms << " ms" << std::endl;
    std::cout << "  Steps: " << result.time.size() << std::endl;
    
    return result;
}

// Write result to CSV
void WriteResultCSV(const SimulationResult& result, const std::string& filename) {
    std::ofstream file(filename);
    file << "time,heave_z,heave_vz\n";
    file << std::setprecision(8);
    for (size_t i = 0; i < result.time.size(); ++i) {
        file << result.time[i] << "," 
             << result.heave_z[i] << "," 
             << result.heave_vz[i] << "\n";
    }
    std::cout << "Wrote: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    // Flush output immediately for debugging
    std::cout.setf(std::ios::unitbuf);
    std::cerr.setf(std::ios::unitbuf);
    
    std::cout << "Sphere Decay Comparison: Convolution vs State-Space\n" << std::endl;
    
    // Initialize Chrono
    SetChronoDataPath(CHRONO_DATA_DIR);
    
    // Paths - try multiple locations
    std::filesystem::path h5_path;
    std::filesystem::path mesh_path;
    
    // Try relative paths from different locations
    std::vector<std::filesystem::path> data_dirs = {
        "../../../../tests/regression/sphere",  // from build/bin/tests/manual
        "../../tests/regression/sphere",
        "../../../tests/regression/sphere",
        "tests/regression/sphere",
        ".",
    };
    
    for (const auto& dir : data_dirs) {
        auto test_h5 = dir / "hydroData" / "sphere.h5";
        auto test_mesh = dir / "geometry" / "oes_task10_sphere.obj";
        if (std::filesystem::exists(test_h5) && std::filesystem::exists(test_mesh)) {
            h5_path = test_h5;
            mesh_path = test_mesh;
            break;
        }
    }
    
    if (h5_path.empty()) {
        std::cerr << "Error: Could not find sphere.h5 and mesh files." << std::endl;
        std::cerr << "Searched in:" << std::endl;
        for (const auto& dir : data_dirs) {
            std::cerr << "  " << std::filesystem::absolute(dir) << std::endl;
        }
        std::cerr << "Please run from the build directory or provide paths." << std::endl;
        return 1;
    }
    
    std::cout << "H5 file: " << std::filesystem::absolute(h5_path) << std::endl;
    std::cout << "Mesh file: " << std::filesystem::absolute(mesh_path) << std::endl;
    
    // Output directory
    std::filesystem::path output_dir = "output";
    ensure_directory_exists(output_dir);
    
    // Run with convolution
    SimulationResult conv_result = RunSphereDecay(
        h5_path.string(),
        mesh_path.string(),
        hydrochrono::hydro::RadiationMethod::kRirfConvolution,
        "RIRF Convolution");
    
    // Run with state-space
    SimulationResult ss_result = RunSphereDecay(
        h5_path.string(),
        mesh_path.string(),
        hydrochrono::hydro::RadiationMethod::kStateSpace,
        "State-Space");
    
    // Write results
    WriteResultCSV(conv_result, (output_dir / "sphere_decay_convolution.csv").string());
    WriteResultCSV(ss_result, (output_dir / "sphere_decay_state_space.csv").string());
    
    // Summary
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Convolution:" << std::endl;
    std::cout << "  Init: " << conv_result.init_time_ms << " ms" << std::endl;
    std::cout << "  Sim:  " << conv_result.sim_time_ms << " ms" << std::endl;
    std::cout << "State-Space:" << std::endl;
    std::cout << "  Init: " << ss_result.init_time_ms << " ms (fitting)" << std::endl;
    std::cout << "  Sim:  " << ss_result.sim_time_ms << " ms" << std::endl;
    
    double sim_speedup = conv_result.sim_time_ms / ss_result.sim_time_ms;
    std::cout << "\nSimulation Loop Speedup: " << sim_speedup << "x" << std::endl;
    
    // Compute max difference in heave position
    double max_diff = 0.0;
    size_t n = std::min(conv_result.heave_z.size(), ss_result.heave_z.size());
    for (size_t i = 0; i < n; ++i) {
        double diff = std::abs(conv_result.heave_z[i] - ss_result.heave_z[i]);
        max_diff = std::max(max_diff, diff);
    }
    std::cout << "Max heave difference: " << max_diff << " m" << std::endl;
    
    // Write kernel fit data
    WriteKernelFit(h5_path.string(), output_dir);
    
    // Write performance scaling data
    WritePerformanceScaling(h5_path.string(), mesh_path.string(), output_dir);
    
    std::cout << "\nRun: python plot_model_comparison.py sphere --save" << std::endl;
    
    return 0;
}

