#include <hydroc/gui/guihelper.h>
#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/assets/ChColor.h>
#include <chrono/core/ChRealtimeStep.h>
#include <chrono>  // std::chrono::high_resolution_clock::now
#include <filesystem>
#include <iomanip>  // std::setprecision
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;

// usage: ./<demos>.exe [DATADIR] [--nogui]
//
// If no argument is given user can set HYDROCHRONO_DATA_DIR
// environment variable to give the data_directory.
//
int main(int argc, char* argv[]) {
    GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";

    if (hydroc::setInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    // Check if --nogui option is set as 2nd argument
    bool visualizationOn = true;
    if (argc > 2 && std::string("--nogui").compare(argv[2]) == 0) {
        visualizationOn = false;
    }

    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfame =
        (DATADIR / "sphere" / "geometry" / "oes_task10_sphere.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "sphere" / "hydroData" / "sphere.h5").lexically_normal().generic_string();

    // system/solver settings
    ChSystemNSC system;
    system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
    double timestep = 0.015;
    system.SetSolverType(ChSolver::Type::GMRES);
    system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
    system.SetStep(timestep);
    ChRealtimeStepTimer realtime_timer;
    double simulationDuration = 600.0;

    // Create user interface
    std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);

    hydroc::gui::UI& ui = *pui.get();

    // Setup Ground
    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetPos(ChVector<>(0, 0, -5));
    ground->SetIdentifier(-1);
    ground->SetBodyFixed(true);
    ground->SetCollide(false);

    // some io/viz options
    bool profilingOn     = true;
    bool saveDataOn      = true;
    std::vector<double> time_vector;
    std::vector<double> heave_position;

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(  //
        body1_meshfame,                                                              // file name
        1000,                                                                        // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // do not collide
    );

    // define the body's initial conditions
    system.Add(sphereBody);
    sphereBody->SetNameString("body1");  // must set body name correctly! (must match .h5 file)
    sphereBody->SetPos(ChVector<>(0, 0, -2));
    sphereBody->SetMass(261.8e3);

    // Create a visualization material
    auto yellow = chrono_types::make_shared<ChVisualMaterial>();
    yellow->SetDiffuseColor(ChColor(0.244f, 0.225f, 0.072f));
    sphereBody->GetVisualShape(0)->SetMaterial(0, yellow);

    // add prismatic joint between sphere and ground (limit to heave motion only)
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(sphereBody, ground, false, ChCoordsys<>(ChVector<>(0, 0, -2)),
                          ChCoordsys<>(ChVector<>(0, 0, -5)));
    system.AddLink(prismatic);

    // Create the spring between body_1 and ground. The spring end points are
    // specified in the body relative frames.
    double rest_length  = 3.0;
    double spring_coef  = 0.0;
    double damping_coef = 0.0;
    auto spring_1       = chrono_types::make_shared<ChLinkTSDA>();
    spring_1->Initialize(sphereBody, ground, false, ChVector<>(0, 0, -2),
                         ChVector<>(0, 0, -5));  // false means positions are in global frame
    // spring_1->SetRestLength(rest_length); // if not set, the rest length is calculated from initial position
    spring_1->SetSpringCoefficient(spring_coef);
    spring_1->SetDampingCoefficient(damping_coef);
    system.AddLink(spring_1);

    auto my_hydro_inputs = std::make_shared<IrregularWave>();
    my_hydro_inputs->wave_height         = 2.0;
    my_hydro_inputs->wave_period         = 12.0;
    my_hydro_inputs->simulation_duration = simulationDuration;
    my_hydro_inputs->simulation_dt       = timestep;
    my_hydro_inputs->ramp_duration       = 60.0;
    // my_hydro_inputs->ramp_duration = 0.0;
    // my_hydro_inputs->SetSpectrumFrequencies(0.001, 1.0, 1000);
    // TODO add option for PiersonMoskowitzSpectrumHz or other spectrum, have a default, do PM for now

    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(sphereBody);
    // TestHydro hydro_forces(bodies, h5fname, my_hydro_inputs);
    TestHydro hydro_forces(bodies, h5fname);
    hydro_forces.AddWaves(my_hydro_inputs);

    // set up free surface from a mesh
    auto fse_plane = chrono_types::make_shared<ChBody>();
    fse_plane->SetPos(ChVector<>(0, 0, 0));
    fse_plane->SetBodyFixed(true);
    fse_plane->SetCollide(false);
    system.AddBody(fse_plane);

    my_hydro_inputs->SetUpWaveMesh();
    std::shared_ptr<ChBody> fse_mesh = chrono_types::make_shared<ChBodyEasyMesh>(  //
        my_hydro_inputs->GetMeshFile(),                                            // file name
        1000,                                                                      // density
        false,                                                                     // do not evaluate mass automatically
        true,                                                                      // create visualization asset
        false                                                                      // do not collide
    );
    fse_mesh->SetMass(1.0);
    fse_mesh->SetPos_dt(my_hydro_inputs->GetWaveMeshVelocity());
    system.Add(fse_mesh);
    auto fse_prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    fse_prismatic->Initialize(fse_plane, fse_mesh, ChCoordsys<>(ChVector<>(1.0, 0.0, 0.0), Q_from_AngY(CH_C_PI_2)));
    system.AddLink(fse_prismatic);

    // Create a visualization material
    auto fse_texture = chrono_types::make_shared<ChVisualMaterial>();
    fse_texture->SetDiffuseColor(ChColor(0.026f, 0.084f, 0.168f));
    fse_texture->SetOpacity(0.1);
    fse_mesh->GetVisualShape(0)->SetMaterial(0, fse_texture);

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();
    // main simulation loop
    ui.Init(&system, "Sphere - Irregular Waves Test");
    ui.SetCamera(8, -25, 15, 0, 0, 0);

    while (system.GetChTime() <= simulationDuration) {
        if (ui.IsRunning(timestep) == false) break;

        if (ui.simulationStarted) {
            system.DoStepDynamics(timestep);

            // append data to output vector
            time_vector.push_back(system.GetChTime());
            heave_position.push_back(sphereBody->GetPos().z());
        }
    }


    // for profiling
    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    if (profilingOn) {
        std::string out_file = "results/sphere_irregular_waves_duration.txt";
        std::ofstream profilingFile(out_file);
        if (!profilingFile.is_open()) {
            if (!std::filesystem::exists("./results")) {
                std::cout << "Path " << std::filesystem::absolute("./results") << " does not exist, creating it now..."
                          << std::endl;
                std::filesystem::create_directories("./results");
                profilingFile.open(out_file);
                if (!profilingFile.is_open()) {
                    std::cout << "Still cannot open file, ending program" << std::endl;
                    return 0;
                }
            }
        }
        profilingFile << duration << "ms\n";
        profilingFile.close();
    }

    if (saveDataOn) {
        std::string out_file = "results/sphere_irreg_waves.txt";
        std::ofstream outputFile(out_file);
        if (!outputFile.is_open()) {
            if (!std::filesystem::exists("./results")) {
                std::cout << "Path " << std::filesystem::absolute("./results") << " does not exist, creating it now..."
                          << std::endl;
                std::filesystem::create_directories("./results");
                outputFile.open(out_file);
                if (!outputFile.is_open()) {
                    std::cout << "Still cannot open file, ending program" << std::endl;
                    return 0;
                }
            }
        }
        outputFile.precision(10);
        outputFile.width(12);
        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(12)
                   << "Heave (m)"
                   //<< std::right << std::setw(18) << "Heave Vel (m/s)"
                   //<< std::right << std::setw(18) << "Heave Force (N)"
                   << std::endl;
        for (int i = 0; i < time_vector.size(); ++i)
            outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
                       << std::right << std::setw(12) << std::setprecision(4) << std::fixed << heave_position[i]
                       << std::endl;
        outputFile.close();
    }

    return 0;
}
