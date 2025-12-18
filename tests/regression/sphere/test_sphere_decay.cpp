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
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    // Initialize environment
    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    // Get model file names
    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfame =
        (DATADIR / "demos" / "sphere" / "geometry" / "oes_task10_sphere.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "demos" / "sphere" / "hydroData" / "sphere.h5").lexically_normal().generic_string();

    // system/solver settings
    ChSystemNSC system;

    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));

    double timestep = 0.015;
    system.SetSolverType(ChSolver::Type::GMRES);
    system.GetSolver()->AsIterative()->SetMaxIterations(300);
    double simulationDuration = 40.0;

    // Output timeseries
    std::vector<double> time_vector;
    std::vector<double> heave_position;

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(
        body1_meshfame,
        1000,   // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // do not collide
    );

    // define the body's initial conditions
    sphereBody->SetName("body1");  // must set body name correctly! (must match .h5 file)
    sphereBody->SetPos(ChVector3d(0, 0, -1));
    sphereBody->SetMass(261.8e3);

    system.Add(sphereBody);

    auto default_dont_add_waves = std::make_shared<NoWave>(1);

    // attach hydrodynamic forces to body
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(sphereBody);

    HydroForces hydro_forces(bodies, h5fname);
    hydro_forces.AddWaves(default_dont_add_waves);

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

    // main simulation loop
    while (system.GetChTime() <= simulationDuration) {
        system.DoStepDynamics(timestep);

        // append data to output vector
        time_vector.push_back(system.GetChTime());
        heave_position.push_back(sphereBody->GetPos().z());
    }

    // for profiling
    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Simulation completed in " << duration << " ms" << std::endl;

    // Save results
    std::string out_dir = hydroc::getTestOutDir() + "/" + RESULTS_DIR_NAME;
    std::filesystem::create_directories(std::filesystem::path(out_dir));

    std::ofstream outputFile(out_dir + "/" + RESULTS_FILE_NAME + ".txt");
    if (outputFile.is_open()) {
        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(12)
                   << "Heave (m)" << std::endl;
        for (size_t i = 0; i < time_vector.size(); ++i)
            outputFile << std::left << std::setw(12) << std::setprecision(6) << std::fixed << time_vector[i]
                       << std::right << std::setw(12) << std::setprecision(6) << std::fixed << heave_position[i]
                       << std::endl;
        outputFile.close();
    } else {
        std::cout << "Error: Could not open output file for writing." << std::endl;
        return 1;
    }

    return 0;
}
