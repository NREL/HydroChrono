#include <hydroc/helper.h>
#include <hydroc/hydro_system.h>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <chrono>
#include <iomanip>
#include <vector>

using namespace chrono;

int main(int argc, char* argv[]) {
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfame =
        (DATADIR / "demos" / "sphere" / "geometry" / "oes_task10_sphere.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "demos" / "sphere" / "hydroData" / "sphere.h5").lexically_normal().generic_string();

    ChSystemNSC system;
    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));

    double timestep = 0.015;
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
    double simulationDuration = hydroc::getSimDuration(40.0, 100.0);

    std::vector<double> time_vector;
    std::vector<double> heave_position;

    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(
        body1_meshfame,
        1000,
        false,
        true,
        false
    );

    sphereBody->SetName("body1");
    sphereBody->SetPos(ChVector3d(0, 0, -1));
    sphereBody->SetMass(261.8e3);
    system.Add(sphereBody);

    auto default_dont_add_waves = std::make_shared<NoWave>();

    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(sphereBody);

    HydroForces hydro_forces(bodies, h5fname);
    hydro_forces.AddWaves(default_dont_add_waves);

    // Use state-space approximation for radiation damping
    hydro_forces.SetRadiationMethod(hydrochrono::hydro::RadiationMethod::kStateSpace);
    hydrochrono::hydro::StateSpaceOptions ss_opts;
    ss_opts.max_order = 10;
    ss_opts.r2_threshold = 0.99;
    hydro_forces.SetStateSpaceOptions(ss_opts);

    auto start = std::chrono::high_resolution_clock::now();

    while (system.GetChTime() <= simulationDuration) {
        system.DoStepDynamics(timestep);
        time_vector.push_back(system.GetChTime());
        heave_position.push_back(sphereBody->GetPos().z());
    }

    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Simulation completed in " << duration << " ms" << std::endl;

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
