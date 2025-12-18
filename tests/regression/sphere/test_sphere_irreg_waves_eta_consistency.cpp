/**
 * @file test_sphere_irreg_waves_eta_consistency.cpp
 * @brief Validates that spectrum-generated and eta-file-imported irregular waves produce identical results.
 *
 * This test ensures the two code paths for irregular wave simulation are consistent:
 * 1. Generate waves from spectrum parameters (wave_height, wave_period)
 * 2. Import waves from an eta file
 *
 * The test runs both simulations with the same parameters and compares the heave response.
 * Any significant difference indicates a bug in one of the code paths.
 */

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

using namespace chrono;

// Shared simulation parameters
const double TIMESTEP      = 0.015;
const double DURATION      = 60.0;  // Short duration for faster testing
const double WAVE_HEIGHT   = 2.0;
const double WAVE_PERIOD   = 12.0;
const double RAMP_DURATION = 0.0;  // No ramp for exact comparison
const int SEED             = 42;   // Fixed seed for reproducibility

int main(int argc, char* argv[]) {
    std::cout << "=== IRREGULAR WAVES ETA CONSISTENCY TEST ===" << std::endl;
    std::cout << "Validates spectrum-generated and eta-imported waves produce identical results.\n" << std::endl;

    // Initialize environment
    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    std::filesystem::path DATADIR(hydroc::getDataDir());
    auto meshfname =
        (DATADIR / "demos" / "sphere" / "geometry" / "oes_task10_sphere.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "demos" / "sphere" / "hydroData" / "sphere.h5").lexically_normal().generic_string();

    // Output directory setup
    std::string out_dir = hydroc::getTestOutDir();
    std::filesystem::create_directories(out_dir + "/" + RESULTS_DIR_NAME);
    std::string eta_file = out_dir + "/" + RESULTS_DIR_NAME + "/temp_eta.txt";

    std::vector<double> heave_spectrum;
    std::vector<double> heave_eta;

    // ========== PHASE 1: Run with spectrum-generated waves ==========
    std::cout << "Phase 1: Running simulation with spectrum-generated waves..." << std::endl;
    {
        ChSystemNSC system;
        system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
        system.SetSolverType(ChSolver::Type::GMRES);
        system.GetSolver()->AsIterative()->SetMaxIterations(300);

        // Ground
        auto ground = chrono_types::make_shared<ChBody>();
        system.AddBody(ground);
        ground->SetPos(ChVector3d(0, 0, -5));
        ground->SetFixed(true);
        ground->EnableCollision(false);

        // Sphere body
        auto sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(meshfname, 1000, false, true, false);
        system.Add(sphereBody);
        sphereBody->SetName("body1");
        sphereBody->SetPos(ChVector3d(0, 0, -2));
        sphereBody->SetMass(261.8e3);

        // Prismatic joint (heave only)
        auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
        prismatic->Initialize(sphereBody, ground, false, ChFramed(ChVector3d(0, 0, -2)),
                              ChFramed(ChVector3d(0, 0, -5)));
        system.AddLink(prismatic);

        // Spring (zero stiffness/damping)
        auto spring = chrono_types::make_shared<ChLinkTSDA>();
        spring->Initialize(sphereBody, ground, false, ChVector3d(0, 0, -2), ChVector3d(0, 0, -5));
        spring->SetSpringCoefficient(0.0);
        spring->SetDampingCoefficient(0.0);
        system.AddLink(spring);

        // Create spectrum-based waves
        IrregularWaveParams spectrum_params;
        spectrum_params.num_bodies_          = 1;
        spectrum_params.simulation_dt_       = TIMESTEP;
        spectrum_params.simulation_duration_ = DURATION;
        spectrum_params.ramp_duration_       = RAMP_DURATION;
        spectrum_params.wave_height_         = WAVE_HEIGHT;
        spectrum_params.wave_period_         = WAVE_PERIOD;
        spectrum_params.frequency_min_       = 0.001;
        spectrum_params.frequency_max_       = 1.0;
        spectrum_params.nfrequencies_        = 1000;
        spectrum_params.seed_                = SEED;

        auto spectrum_waves = std::make_shared<IrregularWaves>(spectrum_params);

        // Setup hydro forces - this initializes the waves via AddH5Data
        std::vector<std::shared_ptr<ChBody>> bodies = {sphereBody};
        HydroForces hydro_forces(bodies, h5fname);
        hydro_forces.AddWaves(spectrum_waves);

        // NOW we can get the free surface data (after AddWaves initialized it)
        auto fse_time      = spectrum_waves->GetFreeSurfaceTime();
        auto fse_elevation = spectrum_waves->GetFreeSurfaceElevation();

        std::cout << "  Generated " << fse_time.size() << " free surface elevation samples." << std::endl;
        if (!fse_time.empty()) {
            std::cout << "  Time range: [" << fse_time.front() << ", " << fse_time.back() << "]" << std::endl;
        }

        // Export eta data to file
        std::cout << "  Exporting eta data to: " << eta_file << std::endl;
        {
            std::ofstream eta_out(eta_file);
            if (!eta_out.is_open()) {
                std::cerr << "ERROR: Could not create eta file: " << eta_file << std::endl;
                return 1;
            }
            eta_out << std::setprecision(10);
            for (size_t i = 0; i < fse_time.size(); ++i) {
                eta_out << fse_time[i] << ":" << fse_elevation[i] << "\n";
            }
            eta_out.close();
        }

        // Run simulation
        while (system.GetChTime() <= DURATION) {
            system.DoStepDynamics(TIMESTEP);
            heave_spectrum.push_back(sphereBody->GetPos().z());
        }
        std::cout << "  Completed. " << heave_spectrum.size() << " timesteps." << std::endl;
    }

    // ========== PHASE 2: Run with eta-file-imported waves ==========
    std::cout << "\nPhase 2: Running simulation with eta-file-imported waves..." << std::endl;
    {
        ChSystemNSC system;
        system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
        system.SetSolverType(ChSolver::Type::GMRES);
        system.GetSolver()->AsIterative()->SetMaxIterations(300);

        // Ground
        auto ground = chrono_types::make_shared<ChBody>();
        system.AddBody(ground);
        ground->SetPos(ChVector3d(0, 0, -5));
        ground->SetFixed(true);
        ground->EnableCollision(false);

        // Sphere body
        auto sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(meshfname, 1000, false, true, false);
        system.Add(sphereBody);
        sphereBody->SetName("body1");
        sphereBody->SetPos(ChVector3d(0, 0, -2));
        sphereBody->SetMass(261.8e3);

        // Prismatic joint (heave only)
        auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
        prismatic->Initialize(sphereBody, ground, false, ChFramed(ChVector3d(0, 0, -2)),
                              ChFramed(ChVector3d(0, 0, -5)));
        system.AddLink(prismatic);

        // Spring (zero stiffness/damping)
        auto spring = chrono_types::make_shared<ChLinkTSDA>();
        spring->Initialize(sphereBody, ground, false, ChVector3d(0, 0, -2), ChVector3d(0, 0, -5));
        spring->SetSpringCoefficient(0.0);
        spring->SetDampingCoefficient(0.0);
        system.AddLink(spring);

        // Create eta-import waves
        IrregularWaveParams eta_params;
        eta_params.num_bodies_          = 1;
        eta_params.simulation_dt_       = TIMESTEP;
        eta_params.simulation_duration_ = DURATION;
        eta_params.ramp_duration_       = RAMP_DURATION;
        eta_params.eta_file_path_       = eta_file;
        eta_params.frequency_min_       = 0.001;
        eta_params.frequency_max_       = 1.0;
        eta_params.nfrequencies_        = 1000;

        auto eta_waves = std::make_shared<IrregularWaves>(eta_params);

        // Setup hydro forces
        std::vector<std::shared_ptr<ChBody>> bodies = {sphereBody};
        HydroForces hydro_forces(bodies, h5fname);
        hydro_forces.AddWaves(eta_waves);

        // Run simulation
        while (system.GetChTime() <= DURATION) {
            system.DoStepDynamics(TIMESTEP);
            heave_eta.push_back(sphereBody->GetPos().z());
        }
        std::cout << "  Completed. " << heave_eta.size() << " timesteps." << std::endl;
    }

    // ========== PHASE 3: Compare results ==========
    std::cout << "\nPhase 3: Comparing results..." << std::endl;

    if (heave_spectrum.size() != heave_eta.size()) {
        std::cerr << "ERROR: Different number of timesteps! Spectrum: " << heave_spectrum.size()
                  << ", Eta: " << heave_eta.size() << std::endl;
        return 1;
    }

    double max_diff     = 0.0;
    double sum_diff_sq  = 0.0;
    size_t max_diff_idx = 0;

    for (size_t i = 0; i < heave_spectrum.size(); ++i) {
        double diff = std::abs(heave_spectrum[i] - heave_eta[i]);
        sum_diff_sq += diff * diff;
        if (diff > max_diff) {
            max_diff     = diff;
            max_diff_idx = i;
        }
    }

    double rms_diff    = std::sqrt(sum_diff_sq / heave_spectrum.size());
    double time_at_max = max_diff_idx * TIMESTEP;

    std::cout << "  Max difference: " << max_diff << " m (at t=" << time_at_max << "s)" << std::endl;
    std::cout << "  RMS difference: " << rms_diff << " m" << std::endl;

    // ========== PHASE 4: Save results ==========
    std::string results_file = out_dir + "/" + RESULTS_DIR_NAME + "/" + RESULTS_FILE_NAME + ".txt";
    std::ofstream out(results_file);
    if (out.is_open()) {
        out << std::setprecision(10);
        out << "Time (s)    Heave_Spectrum (m)    Heave_Eta (m)    Difference (m)\n";
        out << "--------    ------------------    -------------    --------------\n";
        for (size_t i = 0; i < heave_spectrum.size(); ++i) {
            double t = i * TIMESTEP;
            out << std::setw(10) << std::fixed << std::setprecision(3) << t << std::setw(20) << std::fixed
                << std::setprecision(8) << heave_spectrum[i] << std::setw(18) << std::fixed << std::setprecision(8)
                << heave_eta[i] << std::setw(18) << std::scientific << std::setprecision(4)
                << (heave_spectrum[i] - heave_eta[i]) << "\n";
        }
        out.close();
        std::cout << "\n  Results saved to: " << results_file << std::endl;
    }

    // ========== PASS/FAIL determination ==========
    const double tolerance = 1e-6;  // 1 micrometer tolerance
    bool passed            = (max_diff < tolerance);

    std::cout << "\n=== TEST " << (passed ? "PASSED" : "FAILED") << " ===" << std::endl;
    if (!passed) {
        std::cerr << "Spectrum and eta-imported waves produced different results!" << std::endl;
        std::cerr << "Max difference " << max_diff << " exceeds tolerance " << tolerance << std::endl;
        return 1;
    }

    std::cout << "Spectrum-generated and eta-imported waves produce identical results." << std::endl;
    return 0;
}
