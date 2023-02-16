#include <hydroc/hydro_forces.h>

#include <hydroc/helper.h>
//#include "./src/hydro_forces.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono_irrlicht/ChIrrMeshTools.h"
#include "chrono/core/ChRealtimeStep.h"
#include "chrono/physics/ChLinkMate.h"
#include <iomanip>  // std::setprecision
#include <chrono>   // std::chrono::high_resolution_clock::now
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::irrlicht;

// Use the main namespaces of Irrlicht
using namespace irr;
using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;

class MyActionReceiver : public IEventReceiver {
  public:
    MyActionReceiver(ChVisualSystemIrrlicht* vsys, bool& buttonPressed) : pressed(buttonPressed) {
        // store pointer application
        vis = vsys;

        // ..add a GUI button to control pause/play
        pauseButton = vis->GetGUIEnvironment()->addButton(rect<s32>(510, 20, 650, 35));
        buttonText = vis->GetGUIEnvironment()->addStaticText(L"Paused", rect<s32>(560, 20, 600, 35), false);
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

int main(int argc, char* argv[]) {
    // auto start = std::chrono::high_resolution_clock::now();
    GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";

    if (hydroc::setInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfame = (DATADIR / "oswec" / "geometry" / "flap.obj").lexically_normal().generic_string();
    auto body2_meshfame = (DATADIR / "oswec" / "geometry" / "base.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "oswec" / "hydroData" / "oswec.h5").lexically_normal().generic_string();

    // system/solver settings
    ChSystemNSC system;
    system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
    double timestep = 0.03;
    // system.SetTimestepperType(ChTimestepper::Type::HHT);
    system.SetSolverType(ChSolver::Type::GMRES);
    // system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
    system.SetStep(timestep);
    ChRealtimeStepTimer realtime_timer;
    double simulationDuration = 400.0;

    // some io/viz options
    bool visualizationOn = true;
    bool profilingOn = false;
    bool saveDataOn = true;
    std::vector<double> time_vector;
    std::vector<double> flap_rot;

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> flap_body = chrono_types::make_shared<ChBodyEasyMesh>(  //
        body1_meshfame,
        1000,   // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
    std::shared_ptr<ChBody> base_body = chrono_types::make_shared<ChBodyEasyMesh>(  //
        body2_meshfame,
        1000,   // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the float's initial conditions
    system.Add(flap_body);
    flap_body->SetNameString("body1");
    auto ang_rad = CH_C_PI / 18.0;
    flap_body->SetPos(
        ChVector<>(6.1 * std::cos(CH_C_PI / 2.0 - ang_rad), 0.0, -10.0 + 6.1 * std::sin(CH_C_PI / 2.0 - ang_rad)));
    flap_body->SetRot(Q_from_AngAxis(ang_rad, VECT_Y));
    flap_body->SetMass(127000.0);
    flap_body->SetInertiaXX(ChVector<>(1.85e6, 1.85e6, 1.85e6));
    // notes: mass and inertia added to added mass and system mass correctly.

    // define the plate's initial conditions
    system.Add(base_body);
    base_body->SetNameString("body2");
    base_body->SetPos(ChVector<>(0, 0, -10.9));
    base_body->SetMass(999);
    base_body->SetInertiaXX(ChVector<>(1, 1, 1));
    base_body->SetBodyFixed(true);

    // define base-fore flap joint
    ChQuaternion<> revoluteRot = Q_from_AngX(CH_C_PI / 2.0);
    auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
    revolute->Initialize(base_body, flap_body, ChCoordsys<>(ChVector<>(0.0, 0.0, -10.0), revoluteRot));
    system.AddLink(revolute);

    //// define wave parameters (not used in this demo)
    HydroInputs my_hydro_inputs;
    my_hydro_inputs.mode = WaveMode::noWaveCIC;  // or 'regular' or 'regularCIC' or 'irregular';
    ////my_hydro_inputs.regular_wave_amplitude = 0.022;
    ////my_hydro_inputs.regular_wave_omega = 2.10;

    //// attach hydrodynamic forces to body
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(flap_body);
    bodies.push_back(base_body);
    TestHydro blah(bodies, h5fname, my_hydro_inputs);

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

#ifdef HYDROCHRONO_HAVE_IRRLICHT
    if (visualizationOn) {
        // create the irrlicht application for visualizing
        auto irrlichtVis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
        irrlichtVis->AttachSystem(&system);
        irrlichtVis->SetWindowSize(1280, 720);
        irrlichtVis->SetWindowTitle("OSWEC - Decay Test");
        irrlichtVis->SetCameraVertical(CameraVerticalDir::Z);
        irrlichtVis->Initialize();
        irrlichtVis->AddLogo();
        irrlichtVis->AddSkyBox();
        irrlichtVis->AddCamera(ChVector<>(0, -50, -10), ChVector<>(0, 0, -10));  // camera position and where it points
        irrlichtVis->AddTypicalLights();

        // add play/pause button
        bool buttonPressed = false;
        MyActionReceiver receiver(irrlichtVis.get(), buttonPressed);
        irrlichtVis->AddUserEventReceiver(&receiver);

        // main simulation loop
        while (irrlichtVis->Run() && system.GetChTime() <= simulationDuration) {
            irrlichtVis->BeginScene();
            irrlichtVis->Render();
            irrlicht::tools::drawAllCOGs(irrlichtVis.get(), 10.0);
            irrlichtVis->EndScene();
            if (buttonPressed) {
                // step the simulation forwards
                system.DoStepDynamics(timestep);
                // append data to std vector
                time_vector.push_back(system.GetChTime());
                flap_rot.push_back(flap_body->GetRot().Q_to_Euler123().y());
                // force playback to be real-time
                realtime_timer.Spin(timestep);
                // ChSparseMatrix M;
                // system.GetMassMatrix(&M);
                // std::cout << M << std::endl;
            }
        }
    } else {
#endif  // #ifdef HYDROCHRONO_HAVE_IRRLICHT
        int frame = 0;
        while (system.GetChTime() <= simulationDuration) {
            // append data to std vector
            time_vector.push_back(system.GetChTime());
            flap_rot.push_back(flap_body->GetRot().Q_to_Euler123().y());
            // step the simulation forwards
            system.DoStepDynamics(timestep);
            frame++;
        }
#ifdef HYDROCHRONO_HAVE_IRRLICHT
    }
#endif

    // for profiling
    auto end = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    if (profilingOn) {
        std::ofstream profilingFile;
        profilingFile.open("./results/oswec/decay/duration_ms.txt");
        if (!profilingFile.is_open()) {
            if (!std::filesystem::exists("./results/oswec/decay")) {
                std::cout << "Path " << std::filesystem::absolute("./results/oswec/decay")
                          << " does not exist, creating it now..." << std::endl;
                std::filesystem::create_directory("./results");
                std::filesystem::create_directory("./results/oswec");
                std::filesystem::create_directory("./results/oswec/decay");
                profilingFile.open("./results/oswec/decay/duration_ms.txt");
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
        outputFile.open("./results/oswec/decay/oswec_decay.txt");
        if (!outputFile.is_open()) {
            if (!std::filesystem::exists("./results/oswec/decay")) {
                std::cout << "Path " << std::filesystem::absolute("./results/oswec/decay")
                          << " does not exist, creating it now..." << std::endl;
                std::filesystem::create_directory("./results");
                std::filesystem::create_directory("./results/oswec");
                std::filesystem::create_directory("./results/oswec/decay");
                outputFile.open("./results/oswec/decay/oswec_decay.txt");
                if (!outputFile.is_open()) {
                    std::cout << "Still cannot open file, ending program" << std::endl;
                    return 0;
                }
            }
        }
        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(16)
                   << "Flap Rotation y (radians)" << std::right << std::setw(16) << "Flap Rotation y (degrees)"
                   << std::endl;
        for (int i = 0; i < time_vector.size(); ++i)
            outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
                       << std::right << std::setw(16) << std::setprecision(4) << std::fixed << flap_rot[i] << std::right
                       << std::setw(16) << std::setprecision(4) << std::fixed << flap_rot[i] * 360.0 / 6.28
                       << std::endl;
        outputFile.close();
    }
    return 0;
}