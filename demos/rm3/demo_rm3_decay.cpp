#include <hydroc/gui/guihelper.h>
#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/core/ChRealtimeStep.h>

#include <chrono>   // std::chrono::high_resolution_clock::now
#include <iomanip>  // std::setprecision
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;

// usage: ./<demos>.exe [DATADIR] [--nogui]
//
// If no argument is given user can set HYDROCHRONO_DATA_DIR
// environment variable to give the data_directory.
//
int main(int argc, char* argv[]) {
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    SetChronoDataPath(CHRONO_DATA_DIR);

    if (hydroc::SetInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    // Check if --nogui option is set as 2nd argument
    bool visualizationOn = true;
    if (argc > 2 && std::string("--nogui").compare(argv[2]) == 0) {
        visualizationOn = false;
    }

    // Get model file names
    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfame = (DATADIR / "rm3" / "geometry" / "float_cog.obj").lexically_normal().generic_string();
    auto body2_meshfame = (DATADIR / "rm3" / "geometry" / "plate_cog.obj").lexically_normal().generic_string();
    auto h5fname        = (DATADIR / "rm3" / "hydroData" / "rm3.h5").lexically_normal().generic_string();

    // system/solver settings
    ChSystemNSC system;

    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));

    double timestep = 0.01;
    system.SetTimestepperType(ChTimestepper::Type::HHT);
    system.SetSolverType(ChSolver::Type::GMRES);
    system.GetSolver()->AsIterative()->SetMaxIterations(
        300);  // the higher, the easier to keep the constraints satisfied.
    ChRealtimeStepTimer realtime_timer;
    double simulationDuration = 40.0;

    // Create user interface
    std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);

    hydroc::gui::UI& ui = *pui.get();

    // some io/viz options
    bool profilingOn = true;
    bool saveDataOn  = true;
    std::vector<double> time_vector;
    std::vector<double> float_heave_position;
    std::vector<double> plate_heave_position;

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> float_body1 = chrono_types::make_shared<ChBodyEasyMesh>(  //
        body1_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the float's initial conditions
    system.Add(float_body1);
    float_body1->SetName("body1");
    float_body1->SetPos(ChVector3d(0, 0, (-0.72 + 0.1)));
    float_body1->SetMass(725834);
    float_body1->SetInertiaXX(ChVector3d(20907301.0, 21306090.66, 37085481.11));
    // float_body1->EnableCollision(false);

    // Create a visualization material
    auto red = chrono_types::make_shared<ChVisualMaterial>();
    red->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.1f));
    float_body1->GetVisualShape(0)->SetMaterial(0, red);

    // Plate
    std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
    std::shared_ptr<ChBody> plate_body2 = chrono_types::make_shared<ChBodyEasyMesh>(  //
        body2_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // Create a visualization material
    auto blue = chrono_types::make_shared<ChVisualMaterial>();
    blue->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.6f));
    plate_body2->GetVisualShape(0)->SetMaterial(0, blue);

    // define the plate's initial conditions
    system.Add(plate_body2);
    plate_body2->SetName("body2");
    plate_body2->SetPos(ChVector3d(0, 0, (-21.29)));
    plate_body2->SetMass(886691);
    plate_body2->SetInertiaXX(ChVector3d(94419614.57, 94407091.24, 28542224.82));
    // plate_body2->EnableCollision(false);

    // add prismatic joint between the two bodies
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(float_body1, plate_body2, false, ChFramed(ChVector3d(0, 0, -0.72)),
                          ChFramed(ChVector3d(0, 0, -21.29)));
    system.AddLink(prismatic);

    auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
    prismatic_pto->Initialize(float_body1, plate_body2, false, ChVector3d(0, 0, -0.72), ChVector3d(0, 0, -21.29));
    prismatic_pto->SetDampingCoefficient(0.0);
    system.AddLink(prismatic_pto);

    auto default_dont_add_waves = std::make_shared<NoWave>(2);

    // attach hydrodynamic forces to body
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(float_body1);
    bodies.push_back(plate_body2);

    TestHydro hydroForces(bodies, h5fname, default_dont_add_waves);

    //// Debug printing added mass matrix and system mass matrix
    // ChSparseMatrix M;
    // system.GetMassMatrix(&M);
    // std::cout << M << std::endl;

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

    // main simulation loop
    ui.Init(&system, "RM3 - Decay Test");
    ui.SetCamera(0, -50, -10, 0, 0, -10);

    while (system.GetChTime() <= simulationDuration) {
        if (ui.IsRunning(timestep) == false) break;

        if (ui.simulationStarted) {
            system.DoStepDynamics(timestep);

            // append data to output vector
            time_vector.push_back(system.GetChTime());
            float_heave_position.push_back(float_body1->GetPos().z());
            plate_heave_position.push_back(plate_body2->GetPos().z());
        }
    }

    // for profiling
    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    if (profilingOn) {
        std::ofstream profilingFile;
        profilingFile.open("./results/rm3_decay_duration_ms.txt");
        if (!profilingFile.is_open()) {
            if (!std::filesystem::exists("./results")) {
                std::cout << "Path " << std::filesystem::absolute("./results") << " does not exist, creating it now..."
                          << std::endl;
                std::filesystem::create_directory("./results");
                profilingFile.open("./results/rm3_duration_ms.txt");
                if (!profilingFile.is_open()) {
                    std::cout << "Still cannot open file, ending program" << std::endl;
                    return 0;
                }
            }
        }
        profilingFile << duration << "\n";
        profilingFile.close();
    }

    if (saveDataOn) {
        std::ofstream outputFile;
        outputFile.open("./results/rm3_decay.txt");
        if (!outputFile.is_open()) {
            if (!std::filesystem::exists("./results")) {
                std::cout << "Path " << std::filesystem::absolute("./results") << " does not exist, creating it now..."
                          << std::endl;
                std::filesystem::create_directory("./results");
                outputFile.open("./results/rm3_decay.txt");
                if (!outputFile.is_open()) {
                    std::cout << "Still cannot open file, ending program" << std::endl;
                    return 0;
                }
            }
        }
        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(16) << "Float Heave (m)"
                   << std::right << std::setw(16) << "Plate Heave (m)" << std::endl;
        for (int i = 0; i < time_vector.size(); ++i)
            outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
                       << std::right << std::setw(16) << std::setprecision(8) << std::fixed << float_heave_position[i]
                       << std::right << std::setw(16) << std::setprecision(8) << std::fixed << plate_heave_position[i]
                       << std::endl;
        outputFile.close();
    }
    return 0;
}