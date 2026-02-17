#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <chrono>   // std::chrono::high_resolution_clock::now
#include <iomanip>  // std::setprecision
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;

int main(int argc, char* argv[]) {
    std::vector<double> periods = {4.0,  6.0,  8.0,   10.0, 12.0, 14.0, 16.0, 18.0,
                                   18.5, 19.0, 19.25, 19.5, 20.0, 21.0, 22.0, 24.0};
    int reg_wave_num_max        = (int)periods.size();

    // Initialize environment once
    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    // Get model file names
    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfname =
        (DATADIR / "demos" / "oswec" / "geometry" / "flap.obj").lexically_normal().generic_string();
    auto body2_meshfname =
        (DATADIR / "demos" / "oswec" / "geometry" / "base.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "demos" / "oswec" / "hydroData" / "oswec.h5").lexically_normal().generic_string();

    for (int reg_wave_num = 1; reg_wave_num <= reg_wave_num_max; ++reg_wave_num) {
        std::cout << "Wave number: " << reg_wave_num << " of " << reg_wave_num_max << std::endl;

        // system/solver settings
        ChSystemNSC system;

        system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
        double timestep = 0.03;
        system.SetSolverType(ChSolver::Type::GMRES);
        double simulationDuration = hydroc::getSimDuration(1000.0, 1000.0);

        // Output timeseries
        std::vector<double> time_vector;
        std::vector<double> flap_rot;

        // set up body from a mesh
        std::cout << "  Attempting to open mesh file: " << body1_meshfname << std::endl;
        std::shared_ptr<ChBody> flap_body = chrono_types::make_shared<ChBodyEasyMesh>(
            body1_meshfname,
            1000,   // density
            false,  // do not evaluate mass automatically
            true,   // create visualization asset
            false   // collisions
        );

        // define the float's initial conditions
        system.Add(flap_body);
        flap_body->SetName("body1");
        flap_body->SetPos(ChVector3d(0.0, 0.0, -3.9));
        flap_body->SetMass(127000.0);
        flap_body->SetInertiaXX(ChVector3d(1.85e6, 1.85e6, 1.85e6));

        // set up body from a mesh
        std::cout << "  Attempting to open mesh file: " << body2_meshfname << std::endl;
        std::shared_ptr<ChBody> base_body = chrono_types::make_shared<ChBodyEasyMesh>(
            body2_meshfname,
            1000,   // density
            false,  // do not evaluate mass automatically
            true,   // create visualization asset
            false   // collisions
        );

        // define the plate's initial conditions
        system.Add(base_body);
        base_body->SetName("body2");
        base_body->SetPos(ChVector3d(0, 0, -10.15));
        base_body->SetMass(1e9);
        base_body->SetInertiaXX(ChVector3d(1e6, 1e6, 1e6));

        // create ground
        auto ground = chrono_types::make_shared<ChBody>();
        system.AddBody(ground);
        ground->SetPos(ChVector3d(0, 0, -10.15));
        ground->SetTag(-1);
        ground->SetFixed(true);
        ground->EnableCollision(false);

        // fix base to ground with special constraint
        auto anchor = chrono_types::make_shared<ChLinkMateGeneric>();
        anchor->Initialize(base_body, ground, false, base_body->GetVisualModelFrame(),
                           base_body->GetVisualModelFrame());
        system.Add(anchor);
        anchor->SetConstrainedCoords(true, true, true, true, true, true);

        // define base-fore flap joint
        ChQuaternion<> revoluteRot = QuatFromAngleX(CH_PI / 2.0);
        auto revolute              = chrono_types::make_shared<ChLinkLockRevolute>();
        revolute->Initialize(base_body, flap_body, ChFramed(ChVector3d(0.0, 0.0, -8.9), revoluteRot));
        system.AddLink(revolute);

        // attach hydrodynamic forces to body
        std::vector<std::shared_ptr<ChBody>> bodies;
        bodies.push_back(flap_body);
        bodies.push_back(base_body);

        auto my_hydro_inputs = std::make_shared<RegularWave>(static_cast<unsigned int>(bodies.size()));
        my_hydro_inputs->regular_wave_amplitude_ = 0.01;
        my_hydro_inputs->regular_wave_omega_     = (2 * CH_PI) / (periods[reg_wave_num - 1]);

        HydroForces hydro_forces(bodies, h5fname);
        hydro_forces.AddWaves(my_hydro_inputs);

        // for profiling
        auto start = std::chrono::high_resolution_clock::now();

        // main simulation loop
        while (system.GetChTime() < simulationDuration - timestep / 2.0) {
            system.DoStepDynamics(timestep);

            // append data to output vector
            time_vector.push_back(system.GetChTime());
            flap_rot.push_back(flap_body->GetRot().GetCardanAnglesXYZ().y());
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
            outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(12) << "Pitch (rads)"
                       << std::endl;
            for (size_t i = 0; i < time_vector.size(); ++i)
                outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
                           << std::right << std::setw(12) << std::setprecision(4) << std::fixed << flap_rot[i]
                           << std::endl;
            outputFile.close();
        } else {
            std::cout << "Error: Could not open output file for writing." << std::endl;
            return 1;
        }
    }

    return 0;
}
