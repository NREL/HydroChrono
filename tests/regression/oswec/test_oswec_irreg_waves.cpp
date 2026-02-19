#include <chrono>
#include <filesystem>
#include <iomanip>
#include <vector>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <hydroc/helper.h>
#include <hydroc/hydro_system.h>

using namespace chrono;

int main(int argc, char* argv[]) {
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfame = (DATADIR / "demos" / "oswec" / "geometry" / "flap.obj").lexically_normal().generic_string();
    auto body2_meshfame = (DATADIR / "demos" / "oswec" / "geometry" / "base.obj").lexically_normal().generic_string();
    auto h5fname        = (DATADIR / "demos" / "oswec" / "hydroData" / "oswec.h5").lexically_normal().generic_string();

    ChSystemNSC system;
    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
    double timestep = 0.015;
    system.SetSolverType(ChSolver::Type::GMRES);
    double simulationDuration = hydroc::getSimDuration(200.0, 400.0);

    std::vector<double> time_vector;
    std::vector<double> flap_rot;

    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> flap_body = chrono_types::make_shared<ChBodyEasyMesh>(
        body1_meshfame,
        1000,
        false,
        true,
        false
    );

    system.Add(flap_body);
    flap_body->SetName("body1");
    flap_body->SetPos(ChVector3d(0, 0, -3.9));
    flap_body->SetMass(127000.0);
    flap_body->SetInertiaXX(ChVector3d(1.85e6, 1.85e6, 1.85e6));

    std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
    std::shared_ptr<ChBody> base_body = chrono_types::make_shared<ChBodyEasyMesh>(
        body2_meshfame,
        1000,
        false,
        true,
        false
    );

    system.Add(base_body);
    base_body->SetName("body2");
    base_body->SetPos(ChVector3d(0, 0, -10.15));
    base_body->SetMass(999);
    base_body->SetInertiaXX(ChVector3d(1, 1, 1));

    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetPos(ChVector3d(0, 0, -10.15));
    ground->SetTag(-1);
    ground->SetFixed(true);
    ground->EnableCollision(false);

    auto anchor = chrono_types::make_shared<ChLinkMateGeneric>();
    anchor->Initialize(base_body, ground, false, base_body->GetVisualModelFrame(), base_body->GetVisualModelFrame());
    system.Add(anchor);
    anchor->SetConstrainedCoords(true, true, true, true, true, true);

    ChQuaternion<> revoluteRot = QuatFromAngleX(CH_PI / 2.0);
    auto revolute              = chrono_types::make_shared<ChLinkLockRevolute>();
    revolute->Initialize(base_body, flap_body, ChFramed(ChVector3d(0.0, 0.0, -8.9), revoluteRot));
    system.AddLink(revolute);

    // PTO damper between flap and base
    auto pto = chrono_types::make_shared<ChLinkRSDA>();
    pto->Initialize(flap_body, base_body, ChFramed(ChVector3d(0.0, 0.0, -8.9), revoluteRot));
    pto->SetSpringCoefficient(0.0);
    pto->SetDampingCoefficient(12.0e6);
    system.AddLink(pto);

    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(flap_body);
    bodies.push_back(base_body);

    IrregularWaveParams wave_inputs;
    wave_inputs.ramp_duration       = 30.0;
    wave_inputs.wave_height         = 2.5;
    wave_inputs.wave_period         = 8.0;
    wave_inputs.frequency_min       = 0.001;
    wave_inputs.frequency_max       = 1.0;
    wave_inputs.nfrequencies        = 1000;

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

    auto start = std::chrono::high_resolution_clock::now();

    while (system.GetChTime() <= simulationDuration) {
        system.DoStepDynamics(timestep);
        time_vector.push_back(system.GetChTime());
        flap_rot.push_back(flap_body->GetRot().GetCardanAnglesXYZ().y());
    }

    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Simulation completed in " << duration << " ms" << std::endl;

    std::string out_dir = hydroc::getTestOutDir() + "/" + RESULTS_DIR_NAME;
    std::filesystem::create_directories(std::filesystem::path(out_dir));

    std::ofstream outputFile(out_dir + "/" + RESULTS_FILE_NAME + ".txt");
    if (outputFile.is_open()) {
        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(16)
                   << "Flap Rotation y (radians)" << std::right << std::setw(16) << "Flap Rotation y (degrees)"
                   << std::endl;
        for (size_t i = 0; i < time_vector.size(); ++i)
            outputFile << std::left << std::setw(10) << std::setprecision(3) << std::fixed << time_vector[i]
                       << std::right << std::setw(16) << std::setprecision(6) << std::fixed << flap_rot[i]
                       << std::right << std::setw(16) << std::setprecision(6) << std::fixed
                       << flap_rot[i] * 360.0 / (2.0 * CH_PI) << std::endl;
        outputFile.close();
    } else {
        std::cout << "Error: Could not open output file for writing." << std::endl;
        return 1;
    }

    return 0;
}
