/*********************************************************************
 * @file  oswec_decay_ss_comparison.cpp
 * @brief Compare OSWEC decay: Convolution vs State-Space radiation.
 *
 * This diagnostic runs the OSWEC (flap + base) decay test twice:
 *   1. Using RIRF convolution (baseline)
 *   2. Using state-space approximation
 *
 * OSWEC is a 2-body system with 12 DOFs (6 per body), making it a good
 * test case for state-space performance on multi-body systems.
 *
 * Outputs:
 *   - output/oswec_kernel_fit.csv         (t, K_actual, K_fit)
 *   - output/oswec_decay_convolution.csv  (time, flap_pitch_rad)
 *   - output/oswec_decay_state_space.csv  (time, flap_pitch_rad)
 *   - output/oswec_performance.csv        (scaling data)
 *
 * Visualize with: python plot_model_comparison.py oswec --save
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
#include <chrono/physics/ChLinkLock.h>
#include <chrono/physics/ChLinkMate.h>
#include <chrono/solver/ChSolver.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <array>

using namespace chrono;

// =============================================================================
// Physical Constants for OSWEC
// =============================================================================
namespace {
// Geometry
constexpr double kOriginToHingeZ = -8.9;        // m (hinge below origin)
constexpr double kHingeToCgZ = 5.0;             // m (CG above hinge)
constexpr double kBasePositionZ = -10.15;       // m
constexpr double kInitialPitchDeg = 10.0;       // degrees

// Mass properties
constexpr double kFlapMass = 127000.0;          // kg
constexpr double kFlapInertia = 1.85e6;         // kg⋅m²
constexpr double kBaseMass = 999.0;             // kg (nearly fixed)

// Simulation
constexpr double kDefaultTimestep = 0.03;       // s
constexpr double kSimulationDuration = 180.0;   // s (longer for pitch decay)
}  // namespace

// Helper functions for rotation
std::array<double, 3> cross(std::array<double, 3> v1, std::array<double, 3> v2) {
    return {v1[1] * v2[2] - v1[2] * v2[1], 
            v1[2] * v2[0] - v1[0] * v2[2], 
            v1[0] * v2[1] - v1[1] * v2[0]};
}

double dot(std::array<double, 3> v1, std::array<double, 3> v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

std::array<double, 3> normalize(std::array<double, 3> v) {
    double norm = sqrt(dot(v, v));
    return {v[0] / norm, v[1] / norm, v[2] / norm};
}

std::array<double, 3> rotate_vector_3d(std::array<double, 3> vector,
                                       std::array<double, 3> axis,
                                       double angle_in_degrees) {
    double angle_in_radians = angle_in_degrees * M_PI / 180.0;
    axis = normalize(axis);
    std::array<double, 3> rotated_vector;
    for (int i = 0; i < 3; i++) {
        rotated_vector[i] = vector[i] * cos(angle_in_radians) + 
                           cross(axis, vector)[i] * sin(angle_in_radians) +
                           axis[i] * dot(axis, vector) * (1 - cos(angle_in_radians));
    }
    return rotated_vector;
}

std::array<double, 3> add_vectors(std::array<double, 3> v1, std::array<double, 3> v2) {
    return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

void ensure_directory_exists(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path);
    }
}

// Extract and fit flap pitch kernel (body 0, DOF 4 = pitch, self-coupling)
void WriteKernelFit(const std::string& h5_path, const std::filesystem::path& output_dir) {
    std::cout << "\n=== Extracting and Fitting Flap Pitch Kernel ===" << std::endl;
    
    H5FileInfo h5_file(h5_path);
    HydroData hydro_data = h5_file.ReadH5Data();
    
    const int rirf_steps = hydro_data.GetRIRFDims(2);
    const Eigen::VectorXd rirf_time = hydro_data.GetRIRFTimeVector();
    const double dt = (rirf_steps > 1) ? (rirf_time(rirf_steps-1) - rirf_time(0)) / (rirf_steps - 1) : 0.01;
    
    // Flap pitch kernel (body 0, DOF 4 = pitch, column 4)
    const int body = 0;
    const int row_dof = 4;  // pitch
    const int col = 4;      // pitch
    
    Eigen::VectorXd K(rirf_steps);
    for (int step = 0; step < rirf_steps; ++step) {
        K(step) = hydro_data.GetRIRFVal(body, row_dof, col, step);
    }
    
    std::cout << "  RIRF samples: " << rirf_steps << ", dt = " << dt << " s" << std::endl;
    
    hydrochrono::hydro::StateSpaceOptions opts;
    opts.max_order = 15;
    opts.r2_threshold = 0.98;
    opts.max_hankel_size = 400;
    
    hydrochrono::hydro::RadiationStateSpaceFitter fitter(opts);
    hydrochrono::hydro::StateSpaceFitResult result = fitter.FitKernel(K, dt);
    
    std::cout << "  Fit order: " << result.order << ", R² = " << result.r2 << std::endl;
    
    Eigen::VectorXd K_fit = hydrochrono::hydro::RadiationStateSpaceFitter::ReconstructKernel(result, dt, rirf_steps);
    
    std::ofstream file((output_dir / "oswec_kernel_fit.csv").string());
    file << "t,K_actual,K_fit\n" << std::setprecision(8);
    for (int i = 0; i < rirf_steps; ++i) {
        file << rirf_time(i) << "," << K(i) << "," << K_fit(i) << "\n";
    }
    std::cout << "  Wrote: " << (output_dir / "oswec_kernel_fit.csv").string() << std::endl;
}

void WritePerformanceScaling(const std::string& h5_path, const std::string& flap_mesh,
                             const std::string& base_mesh, const std::filesystem::path& output_dir) {
    std::cout << "\n=== Running OSWEC Performance Scaling Test ===" << std::endl;
    
    std::vector<int> step_counts = {100, 500, 1000, 2000};  // Smaller for OSWEC (slower)
    std::vector<double> conv_times, ss_times;
    
    const double dt = 0.03;
    
    // Precompute initial positions
    std::array<double, 3> origin_to_hinge = {0, 0, kOriginToHingeZ};
    std::array<double, 3> hinge_to_cg = {0, 0, kHingeToCgZ};
    std::array<double, 3> axis = {0, 1, 0};
    double angle_in_degrees = kInitialPitchDeg;
    std::array<double, 3> rotated_hinge_to_cg = rotate_vector_3d(hinge_to_cg, axis, angle_in_degrees);
    std::array<double, 3> new_cg = add_vectors(origin_to_hinge, rotated_hinge_to_cg);
    auto ang_rad = CH_PI / 18.0;
    
    for (int steps : step_counts) {
        // Convolution run
        {
            ChSystemNSC system;
            system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
            system.SetSolverType(ChSolver::Type::GMRES);
            
            auto flap = chrono_types::make_shared<ChBodyEasyMesh>(flap_mesh, 1000, false, true, false);
            flap->SetName("body1");
            flap->SetPos(ChVector3d(new_cg[0], new_cg[1], new_cg[2]));
            flap->SetRot(QuatFromAngleY(ang_rad));
            flap->SetMass(kFlapMass);
            flap->SetInertiaXX(ChVector3d(kFlapInertia, kFlapInertia, kFlapInertia));
            system.Add(flap);
            
            auto base = chrono_types::make_shared<ChBodyEasyMesh>(base_mesh, 1000, false, true, false);
            base->SetName("body2");
            base->SetPos(ChVector3d(0, 0, kBasePositionZ));
            base->SetMass(kBaseMass);
            base->SetInertiaXX(ChVector3d(1, 1, 1));
            system.Add(base);
            
            auto ground = chrono_types::make_shared<ChBody>();
            system.AddBody(ground);
            ground->SetPos(ChVector3d(0, 0, kBasePositionZ));
            ground->SetFixed(true);
            
            auto anchor = chrono_types::make_shared<ChLinkMateGeneric>();
            anchor->Initialize(base, ground, false, base->GetVisualModelFrame(), base->GetVisualModelFrame());
            system.Add(anchor);
            anchor->SetConstrainedCoords(true, true, true, true, true, true);
            
            ChQuaternion<> revoluteRot = QuatFromAngleX(CH_PI / 2.0);
            auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
            revolute->Initialize(base, flap, ChFramed(ChVector3d(0.0, 0.0, kOriginToHingeZ), revoluteRot));
            system.AddLink(revolute);
            
            std::vector<std::shared_ptr<ChBody>> bodies = {flap, base};
            HydroSystem hydro(bodies, h5_path);
            hydro.SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kRirfConvolution);
            hydro.AddWaves(std::make_shared<NoWave>(2));
            
            system.DoStepDynamics(dt * 0.001);
            
            auto t1 = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < steps; ++i) system.DoStepDynamics(dt);
            auto t2 = std::chrono::high_resolution_clock::now();
            conv_times.push_back(std::chrono::duration<double, std::micro>(t2 - t1).count());
        }
        
        // State-space run
        {
            ChSystemNSC system;
            system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
            system.SetSolverType(ChSolver::Type::GMRES);
            
            auto flap = chrono_types::make_shared<ChBodyEasyMesh>(flap_mesh, 1000, false, true, false);
            flap->SetName("body1");
            flap->SetPos(ChVector3d(new_cg[0], new_cg[1], new_cg[2]));
            flap->SetRot(QuatFromAngleY(ang_rad));
            flap->SetMass(kFlapMass);
            flap->SetInertiaXX(ChVector3d(kFlapInertia, kFlapInertia, kFlapInertia));
            system.Add(flap);
            
            auto base = chrono_types::make_shared<ChBodyEasyMesh>(base_mesh, 1000, false, true, false);
            base->SetName("body2");
            base->SetPos(ChVector3d(0, 0, kBasePositionZ));
            base->SetMass(kBaseMass);
            base->SetInertiaXX(ChVector3d(1, 1, 1));
            system.Add(base);
            
            auto ground = chrono_types::make_shared<ChBody>();
            system.AddBody(ground);
            ground->SetPos(ChVector3d(0, 0, kBasePositionZ));
            ground->SetFixed(true);
            
            auto anchor = chrono_types::make_shared<ChLinkMateGeneric>();
            anchor->Initialize(base, ground, false, base->GetVisualModelFrame(), base->GetVisualModelFrame());
            system.Add(anchor);
            anchor->SetConstrainedCoords(true, true, true, true, true, true);
            
            ChQuaternion<> revoluteRot = QuatFromAngleX(CH_PI / 2.0);
            auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
            revolute->Initialize(base, flap, ChFramed(ChVector3d(0.0, 0.0, kOriginToHingeZ), revoluteRot));
            system.AddLink(revolute);
            
            std::vector<std::shared_ptr<ChBody>> bodies = {flap, base};
            HydroSystem hydro(bodies, h5_path);
            hydro.SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kStateSpace);
            hydrochrono::hydro::StateSpaceOptions ss_opts;
            ss_opts.max_order = 15;
            ss_opts.r2_threshold = 0.98;
            ss_opts.max_hankel_size = 400;
            hydro.SetStateSpaceOptions(ss_opts);
            hydro.AddWaves(std::make_shared<NoWave>(2));
            
            system.DoStepDynamics(dt * 0.001);  // includes SS fitting
            
            auto t3 = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < steps; ++i) system.DoStepDynamics(dt);
            auto t4 = std::chrono::high_resolution_clock::now();
            ss_times.push_back(std::chrono::duration<double, std::micro>(t4 - t3).count());
        }
        
        std::cout << "  Steps=" << steps << ": Conv=" << conv_times.back()/1000 
                  << "ms, SS=" << ss_times.back()/1000 << "ms" << std::endl;
    }
    
    std::ofstream file((output_dir / "oswec_performance.csv").string());
    file << "sim_steps,conv_time_us,ss_time_us,speedup,conv_per_step_us,ss_per_step_us\n";
    file << std::setprecision(6);
    for (size_t i = 0; i < step_counts.size(); ++i) {
        double speedup = conv_times[i] / ss_times[i];
        double conv_per = conv_times[i] / step_counts[i];
        double ss_per = ss_times[i] / step_counts[i];
        file << step_counts[i] << "," << conv_times[i] << "," << ss_times[i] << ","
             << speedup << "," << conv_per << "," << ss_per << "\n";
    }
    std::cout << "  Wrote: " << (output_dir / "oswec_performance.csv").string() << std::endl;
}

struct SimulationResult {
    std::vector<double> time;
    std::vector<double> flap_pitch_rad;
    std::vector<double> flap_pitch_vel;
    double init_time_ms;
    double sim_time_ms;
    double total_time_ms;
};

SimulationResult RunOswecDecay(
    const std::string& h5_path,
    const std::string& flap_mesh_path,
    const std::string& base_mesh_path,
    hydrochrono::hydro::RadiationMethod method,
    const std::string& method_name) {
    
    std::cout << "\n=== Running OSWEC Decay: " << method_name << " ===" << std::endl;
    
    SimulationResult result;
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // System setup
    ChSystemNSC system;
    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
    
    system.SetSolverType(ChSolver::Type::GMRES);
    
    // Initial rotation geometry
    std::array<double, 3> origin_to_hinge = {0, 0, kOriginToHingeZ};
    std::array<double, 3> hinge_to_cg = {0, 0, kHingeToCgZ};
    std::array<double, 3> axis = {0, 1, 0};
    double angle_in_degrees = kInitialPitchDeg;
    
    std::array<double, 3> rotated_hinge_to_cg = rotate_vector_3d(hinge_to_cg, axis, angle_in_degrees);
    std::array<double, 3> new_cg = add_vectors(origin_to_hinge, rotated_hinge_to_cg);
    
    // Create flap body
    auto flap_body = chrono_types::make_shared<ChBodyEasyMesh>(
        flap_mesh_path, 1000, false, true, false);
    flap_body->SetName("body1");
    auto ang_rad = CH_PI / 18.0;
    flap_body->SetPos(ChVector3d(new_cg[0], new_cg[1], new_cg[2]));
    flap_body->SetRot(QuatFromAngleY(ang_rad));
    flap_body->SetMass(kFlapMass);
    flap_body->SetInertiaXX(ChVector3d(kFlapInertia, kFlapInertia, kFlapInertia));
    system.Add(flap_body);
    
    // Create base body
    auto base_body = chrono_types::make_shared<ChBodyEasyMesh>(
        base_mesh_path, 1000, false, true, false);
    base_body->SetName("body2");
    base_body->SetPos(ChVector3d(0, 0, kBasePositionZ));
    base_body->SetMass(kBaseMass);
    base_body->SetInertiaXX(ChVector3d(1, 1, 1));
    system.Add(base_body);
    
    // Create ground and fix base to it
    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetPos(ChVector3d(0, 0, kBasePositionZ));
    ground->SetTag(-1);
    ground->SetFixed(true);
    ground->EnableCollision(false);
    
    auto anchor = chrono_types::make_shared<ChLinkMateGeneric>();
    anchor->Initialize(base_body, ground, false, base_body->GetVisualModelFrame(), 
                       base_body->GetVisualModelFrame());
    system.Add(anchor);
    anchor->SetConstrainedCoords(true, true, true, true, true, true);
    
    // Revolute joint between base and flap
    ChQuaternion<> revoluteRot = QuatFromAngleX(CH_PI / 2.0);
    auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
    revolute->Initialize(base_body, flap_body, ChFramed(ChVector3d(0.0, 0.0, kOriginToHingeZ), revoluteRot));
    system.AddLink(revolute);
    
    // Setup hydrodynamics
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(flap_body);
    bodies.push_back(base_body);
    
    auto init_start = std::chrono::high_resolution_clock::now();
    
    HydroSystem hydro_forces(bodies, h5_path);
    
    hydro_forces.SetRadiationMethod(method);
    if (method == hydrochrono::hydro::RadiationMethod::kStateSpace) {
        hydrochrono::hydro::StateSpaceOptions opts;
        opts.max_order = 15;           // More modes for complex 2-body system
        opts.r2_threshold = 0.98;      // Higher quality threshold
        opts.max_hankel_size = 400;    // Larger Hankel for better accuracy
        opts.r2_num_samples = 100;     // More samples for R² check
        hydro_forces.SetStateSpaceOptions(opts);
    }
    
    auto no_waves = std::make_shared<NoWave>(2);  // 2 bodies
    hydro_forces.AddWaves(no_waves);
    
    // Force initialization
    system.DoStepDynamics(kDefaultTimestep * 0.001);
    system.SetChTime(0.0);
    flap_body->SetPos(ChVector3d(new_cg[0], new_cg[1], new_cg[2]));
    flap_body->SetRot(QuatFromAngleY(ang_rad));
    flap_body->SetPosDt(ChVector3d(0, 0, 0));
    flap_body->SetRotDt(ChQuaternion<>(1, 0, 0, 0));
    
    auto init_end = std::chrono::high_resolution_clock::now();
    result.init_time_ms = std::chrono::duration<double, std::milli>(init_end - init_start).count();
    
    // Pre-allocate
    int expected_steps = static_cast<int>(kSimulationDuration / kDefaultTimestep) + 10;
    result.time.reserve(expected_steps);
    result.flap_pitch_rad.reserve(expected_steps);
    result.flap_pitch_vel.reserve(expected_steps);
    
    // Simulation loop
    auto sim_start = std::chrono::high_resolution_clock::now();
    
    while (system.GetChTime() <= kSimulationDuration) {
        system.DoStepDynamics(kDefaultTimestep);
        
        result.time.push_back(system.GetChTime());
        result.flap_pitch_rad.push_back(flap_body->GetRot().GetCardanAnglesXYZ().y());
        result.flap_pitch_vel.push_back(flap_body->GetAngVelLocal().y());
    }
    
    auto sim_end = std::chrono::high_resolution_clock::now();
    result.sim_time_ms = std::chrono::duration<double, std::milli>(sim_end - sim_start).count();
    result.total_time_ms = std::chrono::duration<double, std::milli>(sim_end - total_start).count();
    
    std::cout << "  Simulation time: " << kSimulationDuration << " s" << std::endl;
    std::cout << "  Init time: " << result.init_time_ms << " ms" << std::endl;
    std::cout << "  Sim loop time: " << result.sim_time_ms << " ms" << std::endl;
    std::cout << "  Steps: " << result.time.size() << std::endl;
    
    return result;
}

void WriteResultCSV(const SimulationResult& result, const std::string& filename) {
    std::ofstream file(filename);
    file << "time,flap_pitch_rad,flap_pitch_vel\n";
    file << std::setprecision(8);
    for (size_t i = 0; i < result.time.size(); ++i) {
        file << result.time[i] << "," 
             << result.flap_pitch_rad[i] << "," 
             << result.flap_pitch_vel[i] << "\n";
    }
    std::cout << "Wrote: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout.setf(std::ios::unitbuf);
    std::cout << "OSWEC Decay Comparison: Convolution vs State-Space\n" << std::endl;
    
    SetChronoDataPath(CHRONO_DATA_DIR);
    
    // Find data files
    std::filesystem::path h5_path, flap_mesh, base_mesh;
    
    std::vector<std::filesystem::path> data_dirs = {
        "../../../../tests/regression/oswec",
        "../../tests/regression/oswec",
        "../../../tests/regression/oswec",
        "tests/regression/oswec",
    };
    
    for (const auto& dir : data_dirs) {
        auto test_h5 = dir / "hydroData" / "oswec.h5";
        auto test_flap = dir / "geometry" / "flap.obj";
        auto test_base = dir / "geometry" / "base.obj";
        if (std::filesystem::exists(test_h5) && 
            std::filesystem::exists(test_flap) && 
            std::filesystem::exists(test_base)) {
            h5_path = test_h5;
            flap_mesh = test_flap;
            base_mesh = test_base;
            break;
        }
    }
    
    if (h5_path.empty()) {
        std::cerr << "Error: Could not find OSWEC data files." << std::endl;
        return 1;
    }
    
    std::cout << "H5 file: " << std::filesystem::absolute(h5_path) << std::endl;
    std::cout << "Flap mesh: " << std::filesystem::absolute(flap_mesh) << std::endl;
    std::cout << "Base mesh: " << std::filesystem::absolute(base_mesh) << std::endl;
    
    std::filesystem::path output_dir = "output";
    ensure_directory_exists(output_dir);
    
    // Run convolution
    SimulationResult conv_result = RunOswecDecay(
        h5_path.string(), flap_mesh.string(), base_mesh.string(),
        hydrochrono::hydro::RadiationMethod::kRirfConvolution,
        "RIRF Convolution");
    
    // Run state-space
    SimulationResult ss_result = RunOswecDecay(
        h5_path.string(), flap_mesh.string(), base_mesh.string(),
        hydrochrono::hydro::RadiationMethod::kStateSpace,
        "State-Space");
    
    // Write results
    WriteResultCSV(conv_result, (output_dir / "oswec_decay_convolution.csv").string());
    WriteResultCSV(ss_result, (output_dir / "oswec_decay_state_space.csv").string());
    
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
    
    // Max difference
    double max_diff = 0.0;
    size_t n = std::min(conv_result.flap_pitch_rad.size(), ss_result.flap_pitch_rad.size());
    for (size_t i = 0; i < n; ++i) {
        double diff = std::abs(conv_result.flap_pitch_rad[i] - ss_result.flap_pitch_rad[i]);
        max_diff = std::max(max_diff, diff);
    }
    std::cout << "Max pitch difference: " << (max_diff * 180.0 / M_PI) << " degrees" << std::endl;
    
    // Write kernel fit data
    WriteKernelFit(h5_path.string(), output_dir);
    
    // Write performance scaling data
    WritePerformanceScaling(h5_path.string(), flap_mesh.string(), base_mesh.string(), output_dir);
    
    std::cout << "\nRun: python plot_model_comparison.py oswec --save" << std::endl;
    
    return 0;
}

