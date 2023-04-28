#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#ifdef HYDROCHRONO_HAVE_IRRLICHT
    #include "chrono_irrlicht/ChIrrMeshTools.h"
    #include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
// Use the main namespaces of Irrlicht
using namespace irr;
using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;
using namespace chrono::irrlicht;
#endif

#include <chrono/core/ChRealtimeStep.h>
#include <chrono/assets/ChColor.h>
#include <chrono>  // std::chrono::high_resolution_clock::now
#include <filesystem>
#include <iomanip>  // std::setprecision
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;

#ifdef HYDROCHRONO_HAVE_IRRLICHT
class MyActionReceiver : public IEventReceiver {
  public:
    MyActionReceiver(ChVisualSystemIrrlicht* vsys, bool& buttonPressed) : pressed(buttonPressed) {
        // store pointer application
        vis = vsys;

        // ..add a GUI button to control pause/play
        pauseButton = vis->GetGUIEnvironment()->addButton(rect<s32>(510, 20, 650, 35));
        buttonText  = vis->GetGUIEnvironment()->addStaticText(L"Paused", rect<s32>(560, 20, 600, 35), false);
    }

    bool OnEvent(const SEvent& event) {
        // check if user clicked button
        if (event.EventType == EET_GUI_EVENT) {
            switch (event.GUIEvent.EventType) {
                case EGET_BUTTON_CLICKED:
                    pressed = !pressed;
                    if (pressed) {
                        buttonText->setText(L"Playing");
                    } else {
                        buttonText->setText(L"Paused");
                    }
                    return pressed;
                    break;
                default:
                    break;
            }
        }
        return false;
    }

  private:
    ChVisualSystemIrrlicht* vis;
    IGUIButton* pauseButton;
    IGUIStaticText* buttonText;

    bool& pressed;
};
#endif

int main(int argc, char* argv[]) {
    GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";

    if (hydroc::setInitialEnvironment(argc, argv) != 0) {
        return 1;
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

    // Setup Ground
    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetPos(ChVector<>(0, 0, -5));
    ground->SetIdentifier(-1);
    ground->SetBodyFixed(true);
    ground->SetCollide(false);

    // some io/viz options
    bool visualizationOn = true;
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
    //my_hydro_inputs->mode                   = WaveMode::irregular;                     // uses regular wave mode
    my_hydro_inputs->wave_height            = 2.0;
    my_hydro_inputs->wave_period            = 12.0;
    my_hydro_inputs->simulation_duration    = simulationDuration;
    my_hydro_inputs->simulation_dt          = timestep;
    my_hydro_inputs->ramp_duration          = 60.0;
    //my_hydro_inputs->ramp_duration = 0.0;
    //my_hydro_inputs->SetSpectrumFrequencies(0.001, 1.0, 1000);
    //TODO add option for PiersonMoskowitzSpectrumHz or other spectrum, have a default, do PM for now

    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(sphereBody);
    //TestHydro hydro_forces(bodies, h5fname, my_hydro_inputs);
    TestHydro hydro_forces(bodies, h5fname);
    hydro_forces.AddWaves(my_hydro_inputs);

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

    // set up free surface from a mesh 
    auto fse_plane = chrono_types::make_shared<ChBody>();
    fse_plane->SetPos(ChVector<>(0, 0, 0));
    fse_plane->SetBodyFixed(true);
    fse_plane->SetCollide(false);
    system.AddBody(fse_plane);

    my_hydro_inputs->SetUpWaveMesh();
    std::shared_ptr<ChBody> fse_mesh = chrono_types::make_shared<ChBodyEasyMesh>(  //
        my_hydro_inputs->GetMeshFile(),                                                              // file name
        1000,                                                                        // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // do not collide
    );
    fse_mesh->SetMass(1.0);
    fse_mesh->SetPos_dt(my_hydro_inputs->GetWaveMeshVelocity());
    system.Add(fse_mesh);
    auto fse_prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    fse_prismatic->Initialize(fse_plane, fse_mesh, ChCoordsys<>(ChVector<>(1.0, 0.0, 0.0), Q_from_AngY(CH_C_PI_2)));
    system.AddLink(fse_prismatic);

#ifdef HYDROCHRONO_HAVE_IRRLICHT
    if (visualizationOn) {
        // Create a visualization material
        auto yellow = chrono_types::make_shared<ChVisualMaterial>();
        yellow->SetDiffuseColor(ChColor(0.244f, 0.225f, 0.072f));
        sphereBody->GetVisualShape(0)->SetMaterial(0, yellow);

        // Create a visualization material
        auto fse_texture = chrono_types::make_shared<ChVisualMaterial>();
        fse_texture->SetDiffuseColor(ChColor(0.026f, 0.084f, 0.168f));
        fse_texture->SetOpacity(0.1);
        fse_mesh->GetVisualShape(0)->SetMaterial(0, fse_texture);

        // create the irrlicht application for visualizing
        auto irrlichtVis = chrono_types::make_shared<ChVisualSystemIrrlicht>();

        irrlichtVis->AttachSystem(&system);
        irrlichtVis->SetWindowSize(1280, 720);
        irrlichtVis->SetWindowTitle("Sphere - Irregular Waves Test");
        irrlichtVis->SetCameraVertical(CameraVerticalDir::Z);
        irrlichtVis->Initialize();

        irrlichtVis->AddLogo();
        irrlichtVis->AddSkyBox();
        irrlichtVis->AddCamera(ChVector<>(8, -25, 15), ChVector<>(0, 0, 0));
        irrlichtVis->AddTypicalLights();

        // add play/pause button
        bool buttonPressed = false;
        MyActionReceiver receiver(irrlichtVis.get(), buttonPressed);
        irrlichtVis->AddUserEventReceiver(&receiver);

        // Main simulation loop
        int frame = 0;
        while (irrlichtVis->Run() && system.GetChTime() <= simulationDuration) {
            irrlichtVis->BeginScene();
            irrlichtVis->Render();

            // Add grid to materialize horizontal plane
            tools::drawGrid(irrlichtVis.get(), 1, 1, 30, 30,
                            ChCoordsys<>(ChVector<>(0, 0.0, 0), Q_from_AngZ(CH_C_PI_2)), chrono::ChColor(.1f, .1f, .1f),
                            true);

            if (buttonPressed) {
                // step simulation forward
                system.DoStepDynamics(timestep);
                // append data to std vector
                time_vector.push_back(system.GetChTime());
                heave_position.push_back(sphereBody->GetPos().z());
                frame++;
            }
            irrlichtVis->EndScene();
        }
    } else {
#endif  // #ifdef HYDROCHRONO_HAVE_IRRLICHT
        int frame = 0;
        while (system.GetChTime() <= simulationDuration) {
            // step the simulation forwards
            system.DoStepDynamics(timestep);
            // append data to std vector
            time_vector.push_back(system.GetChTime());
            heave_position.push_back(sphereBody->GetPos().z());
            frame++;
        }
#ifdef HYDROCHRONO_HAVE_IRRLICHT
    }
#endif

    // for profiling
    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    if (profilingOn) {
        std::string out_dir  = "results/sphere_irregular_waves/";
        std::string out_file = "duration_ms.txt";
        std::ofstream profilingFile(out_dir + out_file);
        if (!profilingFile.is_open()) {
            if (!std::filesystem::exists(out_dir)) {
                std::cout << "Path " << std::filesystem::absolute(out_dir) << " does not exist, creating it now..."
                          << std::endl;
                std::filesystem::create_directories(out_dir);
                profilingFile.open(out_dir + out_file);
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
        std::string out_dir  = "results/sphere_irregular_waves/";
        std::string out_file = "irreg_test.txt"; //"irreg_H_" + my_hydro_inputs.wave_height + "_T_" + my_hydro_inputs.wave_period + ".txt ";
        std::ofstream outputFile(out_dir + out_file);
        if (!outputFile.is_open()) {
            if (!std::filesystem::exists(out_dir)) {
                std::cout << "Path " << std::filesystem::absolute(out_dir) << " does not exist, creating it now..."
                          << std::endl;
                std::filesystem::create_directories(out_dir);
                outputFile.open(out_dir + out_file);
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
