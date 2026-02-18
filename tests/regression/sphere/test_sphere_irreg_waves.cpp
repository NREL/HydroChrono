#include <chrono>
#include <filesystem>
#include <iomanip>
#include <vector>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

// Use the namespaces of Chrono
using namespace chrono;

int main(int argc, char* argv[]) {
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    // Initialize environment
    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfame =
        (DATADIR / "demos" / "sphere" / "geometry" / "oes_task10_sphere.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "demos" / "sphere" / "hydroData" / "sphere.h5").lexically_normal().generic_string();

    // system/solver settings
    ChSystemNSC system;
    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
    double timestep = 0.015;
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
    double simulationDuration = hydroc::getSimDuration(600.0, 1200.0);

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
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(
        body1_meshfame,
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

    std::cout << "Body created from the mesh file: " << body1_meshfame << std::endl;

    // add prismatic joint between sphere and ground (limit to heave motion only)
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(sphereBody, ground, false, ChFramed(ChVector3d(0, 0, -2)), ChFramed(ChVector3d(0, 0, -5)));
    system.AddLink(prismatic);

    // create the spring between body_1 and ground
    double spring_coef  = 0.0;
    double damping_coef = 0.0;
    auto spring_1       = chrono_types::make_shared<ChLinkTSDA>();
    spring_1->Initialize(sphereBody, ground, false, ChVector3d(0, 0, -2), ChVector3d(0, 0, -5));
    spring_1->SetSpringCoefficient(spring_coef);
    spring_1->SetDampingCoefficient(damping_coef);
    system.AddLink(spring_1);

    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(sphereBody);

    IrregularWaveParams wave_inputs;
    wave_inputs.simulation_dt_       = timestep;
    wave_inputs.simulation_duration_ = simulationDuration;
    wave_inputs.ramp_duration_       = 60.0;
    wave_inputs.wave_height_         = 2.0;
    wave_inputs.wave_period_         = 12.0;
    wave_inputs.frequency_min_       = 0.001;
    wave_inputs.frequency_max_       = 1.0;
    wave_inputs.nfrequencies_        = 1000;

    std::shared_ptr<IrregularWaves> my_hydro_inputs;

    try {
        my_hydro_inputs = std::make_shared<IrregularWaves>(wave_inputs);
    } catch (const std::exception& e) {
        std::cerr << "Caught exception: " << e.what() << '\n';
        return 1;
    } catch (...) {
        std::cerr << "Caught unknown exception.\n";
        return 1;
    }

    if (!my_hydro_inputs) {
        std::cerr << "ERROR: Failed to create IrregularWaves object." << std::endl;
        return 1;
    }

    HydroForces hydro_forces(bodies, h5fname);
    hydro_forces.AddWaves(my_hydro_inputs);

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
        outputFile.precision(10);
        outputFile.width(12);
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

    return 0;
}
