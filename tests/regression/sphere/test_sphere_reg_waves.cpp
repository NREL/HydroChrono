#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <chrono>  // std::chrono::high_resolution_clock::now
#include <filesystem>
#include <iomanip>  // std::setprecision
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;

int main(int argc, char* argv[]) {
    std::vector<double> task10_wave_amps_0005 = {0.044, 0.078, 0.095, 0.123, 0.177, 0.24, 0.314, 0.397, 0.491, 0.594};
    std::vector<double> task10_wave_amps_002  = {0.177, 0.314, 0.380, 0.491, 0.706, 0.961, 1.256, 1.589, 1.962, 2.374};
    std::vector<double> task10_wave_amps      = task10_wave_amps_002;

    double task10_wave_omegas[]    = {2.094395102, 1.570796327, 1.427996661, 1.256637061, 1.047197551,
                                      0.897597901, 0.785398163, 0.698131701, 0.628318531, 0.571198664};
    double task10_damping_coeffs[] = {398736.034, 118149.758, 90080.857,  161048.558, 322292.419,
                                      479668.979, 633979.761, 784083.286, 932117.647, 1077123.445};
    int reg_wave_num_max           = (int)task10_wave_amps.size();

    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    // Initialize environment
    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    // Get model file names
    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfname =
        (DATADIR / "demos" / "sphere" / "geometry" / "oes_task10_sphere.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "demos" / "sphere" / "hydroData" / "sphere.h5").lexically_normal().generic_string();

    for (int reg_wave_num = 1; reg_wave_num <= reg_wave_num_max; ++reg_wave_num) {
        std::cout << "Wave number: " << reg_wave_num << " of " << reg_wave_num_max << std::endl;

        // system/solver settings
        ChSystemNSC system;
        system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
        double timestep = 0.015;
        system.SetSolverType(ChSolver::Type::SPARSE_QR);
        double simulation_duration = hydroc::getSimDuration(240.0, 800.0);

        // Setup Ground
        auto ground = chrono_types::make_shared<ChBody>();
        system.AddBody(ground);
        ground->SetPos(ChVector3d(0, 0, -5));
        ground->SetTag(-1);
        ground->SetFixed(true);
        ground->EnableCollision(false);

        // Output timeseries
        std::vector<double> time_vector;
        std::vector<double> heave_position;

        // set up body from a mesh
        std::cout << "  Attempting to open mesh file: " << body1_meshfname << std::endl;
        std::shared_ptr<ChBody> sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(
            body1_meshfname,
            1000,   // density
            false,  // do not evaluate mass automatically
            true,   // create visualization asset
            false   // do not collide
        );

        // define the body's initial conditions
        system.Add(sphereBody);
        sphereBody->SetName("body1");  // must set body name correctly! (must match .h5 file)
        sphereBody->SetPos(ChVector3d(0, 0, -2));
        sphereBody->SetMass(261.8e3);

        // add prismatic joint between sphere and ground (limit to heave motion only)
        auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
        prismatic->Initialize(sphereBody, ground, false, ChFramed(ChVector3d(0, 0, -2)),
                              ChFramed(ChVector3d(0, 0, -5)));
        system.AddLink(prismatic);

        // Create the spring between body_1 and ground
        double spring_coef  = 0.0;
        double damping_coef = task10_damping_coeffs[reg_wave_num - 1];
        auto spring_1       = chrono_types::make_shared<ChLinkTSDA>();
        spring_1->Initialize(sphereBody, ground, false, ChVector3d(0, 0, -2), ChVector3d(0, 0, -5));
        spring_1->SetSpringCoefficient(spring_coef);
        spring_1->SetDampingCoefficient(damping_coef);
        system.AddLink(spring_1);

        auto my_hydro_inputs                     = std::make_shared<RegularWave>(1);
        my_hydro_inputs->regular_wave_amplitude_ = task10_wave_amps[reg_wave_num - 1];
        my_hydro_inputs->regular_wave_omega_     = task10_wave_omegas[reg_wave_num - 1];

        std::vector<std::shared_ptr<ChBody>> bodies;
        bodies.push_back(sphereBody);
        HydroForces hydro_forces(bodies, h5fname);
        hydro_forces.AddWaves(my_hydro_inputs);

        // for profiling
        auto start = std::chrono::high_resolution_clock::now();

        // main simulation loop
        while (system.GetChTime() < simulation_duration - timestep / 2.0) {
            system.DoStepDynamics(timestep);

            // append data to output vector
            time_vector.push_back(system.GetChTime());
            heave_position.push_back(sphereBody->GetPos().z());
        }

        // for profiling
        auto end          = std::chrono::high_resolution_clock::now();
        unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "  Simulation completed in " << duration << " ms" << std::endl;

        // Save results
        std::string out_dir = hydroc::getTestOutDir() + "/" + RESULTS_DIR_NAME;
        std::filesystem::create_directories(std::filesystem::path(out_dir));

        std::ofstream outputFile(out_dir + "/" + RESULTS_FILE_NAME + "_" + std::to_string(reg_wave_num) + ".txt");
        if (outputFile.is_open()) {
            outputFile.precision(10);
            outputFile.width(12);
            outputFile << "Wave #: \t" << reg_wave_num << "\n";
            outputFile << "Wave amplitude (m): \t" << my_hydro_inputs->regular_wave_amplitude_ << "\n";
            outputFile << "Wave omega (rad/s): \t" << my_hydro_inputs->regular_wave_omega_ << "\n";
            outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(12) << "Heave (m)" << "\n";
            outputFile << std::left << std::setw(10) << "----------" << std::right << std::setw(12) << "----------" << "\n";

            for (size_t i = 0; i < time_vector.size(); i++) {
                outputFile << std::left << std::setw(10) << std::fixed << std::setprecision(3) << time_vector[i]
                           << std::right << std::setw(12) << std::fixed << std::setprecision(6) << heave_position[i]
                           << "\n";
            }
            outputFile.close();
        } else {
            std::cout << "Error: Could not open output file for writing." << std::endl;
            return 1;
        }

        // Clear vectors for next iteration
        time_vector.clear();
        heave_position.clear();
    }

    return 0;
}
