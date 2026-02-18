#include <hydroc/helper.h>
#include <hydroc/hydro_system.h>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemSMC.h>

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

    auto body1_meshfame = (DATADIR / "demos" / "f3of" / "geometry" / "base.obj").lexically_normal().generic_string();
    auto body2_meshfame = (DATADIR / "demos" / "f3of" / "geometry" / "flap.obj").lexically_normal().generic_string();
    auto body3_meshfame = (DATADIR / "demos" / "f3of" / "geometry" / "flap.obj").lexically_normal().generic_string();
    auto h5fname        = (DATADIR / "demos" / "f3of" / "hydroData" / "f3of.h5").lexically_normal().generic_string();

    // system/solver settings
    ChSystemSMC system;

    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
    double timestep = 0.02;
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
    double simulationDuration = hydroc::getSimDuration(100.0, 240.0);

    // Output timeseries
    std::vector<double> time_vector;
    std::vector<double> base_surge;
    std::vector<double> base_pitch;
    std::vector<double> fore_pitch;
    std::vector<double> aft_pitch;

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> base = chrono_types::make_shared<ChBodyEasyMesh>(
        body1_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the base's initial conditions
    system.Add(base);
    base->SetName("body1");
    base->SetMass(1089825.0);
    base->SetInertiaXX(ChVector3d(100000000.0, 76300000.0, 100000000.0));

    std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
    std::shared_ptr<ChBody> flapFore = chrono_types::make_shared<ChBodyEasyMesh>(
        body2_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the fore flap's initial conditions
    system.Add(flapFore);
    flapFore->SetName("body2");
    flapFore->SetMass(179250.0);
    flapFore->SetInertiaXX(ChVector3d(100000000.0, 1300000.0, 100000000.0));

    std::cout << "Attempting to open mesh file: " << body3_meshfame << std::endl;
    std::shared_ptr<ChBody> flapAft = chrono_types::make_shared<ChBodyEasyMesh>(
        body3_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the aft flap's initial conditions
    system.Add(flapAft);
    flapAft->SetName("body3");
    flapAft->SetMass(179250.0);
    flapAft->SetInertiaXX(ChVector3d(100000000.0, 1300000.0, 100000000.0));

    // DT1 set up (surge decay, flaps locked, no waves)
    base->SetPos(ChVector3d(5.0, 0.0, -9.0));
    flapFore->SetPos(ChVector3d(5.0 + -12.5, 0.0, -9.0 + 3.5));
    flapAft->SetPos(ChVector3d(5.0 + 12.5, 0.0, -9.0 + 3.5));

    // set up revolute joints and lock them
    auto revoluteFore          = chrono_types::make_shared<ChLinkLockRevolute>();
    auto revoluteAft           = chrono_types::make_shared<ChLinkLockRevolute>();
    ChQuaternion<> revoluteRot = QuatFromAngleX(CH_PI / 2.0);
    revoluteFore->Initialize(base, flapFore, ChFramed(ChVector3d(5.0 - 12.5, 0.0, -9.0), revoluteRot));
    system.AddLink(revoluteFore);
    revoluteAft->Initialize(base, flapAft, ChFramed(ChVector3d(5.0 + 12.5, 0.0, -9.0), revoluteRot));
    system.AddLink(revoluteAft);
    revoluteFore->Lock(true);
    revoluteAft->Lock(true);

    // create ground
    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetPos(ChVector3d(0, 0, -9.0));
    ground->SetTag(-1);
    ground->SetFixed(true);
    ground->EnableCollision(false);

    // add prismatic joint between the base and ground
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(ground, base, ChFramed(ChVector3d(0.0, 0.0, -9.0), QuatFromAngleY(CH_PI_2)));
    system.AddLink(prismatic);

    // add damping to prismatic joint
    auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
    prismatic_pto->Initialize(ground, base, true, ChVector3d(0.0, 0.0, 0.0), ChVector3d(0.0, 0.0, 0.0));
    prismatic_pto->SetSpringCoefficient(1e5);
    prismatic_pto->SetRestLength(0.0);
    system.AddLink(prismatic_pto);

    // define wave parameters
    auto default_dont_add_waves = std::make_shared<NoWave>();

    // set up hydro forces
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(base);
    bodies.push_back(flapFore);
    bodies.push_back(flapAft);
    HydroForces hydroforces(bodies, h5fname, default_dont_add_waves);

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

    // main simulation loop
    while (system.GetChTime() <= simulationDuration) {
        system.DoStepDynamics(timestep);

        // append data to output vector
        time_vector.push_back(system.GetChTime());
        base_surge.push_back(base->GetPos().x());
        base_pitch.push_back(base->GetRot().GetCardanAnglesXYZ().y());
        fore_pitch.push_back(flapFore->GetRot().GetCardanAnglesXYZ().y());
        aft_pitch.push_back(flapAft->GetRot().GetCardanAnglesXYZ().y());
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
        outputFile
            << "Time (s)    Base Surge (m)Base Pitch (radians)Flap Fore Pitch (radians)Flap Aft Pitch (radians)"
            << std::endl;
        for (size_t i = 0; i < time_vector.size(); i++) {
            outputFile << std::fixed << std::setprecision(4) << time_vector[i] << "                "
                       << base_surge[i] << "         " << base_pitch[i] << "         " << fore_pitch[i]
                       << "          " << aft_pitch[i] << std::endl;
        }
        outputFile.close();
    } else {
        std::cout << "Error: Could not open output file for writing." << std::endl;
        return 1;
    }

    return 0;
}
