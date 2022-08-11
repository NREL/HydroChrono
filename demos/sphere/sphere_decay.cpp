//#include "../../src/hydro_forces.h"
#include "./src/hydro_forces.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono_irrlicht/ChIrrMeshTools.h"
#include "chrono/core/ChRealtimeStep.h"
#include <iomanip> // std::setprecision

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

// Define a class to manage user inputs via the GUI (i.e. play/pause button)

class MyActionReceiver : public IEventReceiver {
	public:
		MyActionReceiver(ChVisualSystemIrrlicht* vsys, bool& buttonPressed)
			: pressed(buttonPressed) {
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
					}
					else {
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

// the main program to be executed:

int main(int argc, char* argv[]) {
	//auto start = std::chrono::high_resolution_clock::now();
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	// system/solver settings
	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0, 0, -9.81));
	double timestep = 0.06;
	system.SetSolverType(ChSolver::Type::GMRES);
	system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
	system.SetStep(timestep);
	ChRealtimeStepTimer realtime_timer;
	bool visualizationOn = true;
	double simulationDuration = 40.0;
	//double simulationStartTime = 0.0;

	// set up body from a mesh
	std::shared_ptr<ChBody> sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(       //
		GetChronoDataFile("../../HydroChrono/demos/sphere/geometry/oes_task10_sphere.obj").c_str(),   // file name
		1000,                                                                             // density
		false,                                                                            // do not evaluate mass automatically
		true,                                                                             // create visualization asset
		false,                                                                            // do not collide
		nullptr,                                                                          // no need for contact material
		0                                                                                 // swept sphere radius
		);

	// define the body's initial conditions
	system.Add(sphereBody);
	sphereBody->SetNameString("body1"); // must set body name correctly! (must match .h5 file)
	sphereBody->SetPos(ChVector<>(0, 0, -1));
	sphereBody->SetMass(261.8e3);

	// define wave parameters (not used in this demo)
	HydroInputs my_hydro_inputs;
	my_hydro_inputs.SetRegularWaveAmplitude(0.022);
	my_hydro_inputs.SetRegularWaveOmega(2.10);

	// attach hydrodynamic forces to body
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(sphereBody);
	TestHydro blah(bodies, "../../HydroChrono/demos/sphere/hydroData/sphere.h5", my_hydro_inputs);

	// set up output file
	std::string outputFileName = "sphere_decay.txt";                    /// < put name of your output file here
	std::ofstream zPosition(outputFileName, std::ofstream::out);
	if (!zPosition.is_open()) {
		std::cout << "Error opening file \"" + outputFileName + "\". Please make sure this file path exists then try again\n";
		return -1;
	}

	zPosition << std::left << std::setw(10) << "Time (s)"
		<< std::right << std::setw(12) << "Heave(m)"
		<< std::right << std::setw(18) << "Heave Vel (m/s)" 
		<< std::right << std::setw(18) << "Heave Force(N)"
		<< std::endl;

	if (visualizationOn){
		// create the irrlicht application for visualizing
		auto irrlichtVis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
		irrlichtVis->AttachSystem(&system);
		irrlichtVis->SetWindowSize(1280, 720);
		irrlichtVis->SetWindowTitle("Sphere - Decay Test");
		irrlichtVis->SetCameraVertical(CameraVerticalDir::Z);
		irrlichtVis->Initialize();
		irrlichtVis->AddLogo();
		irrlichtVis->AddSkyBox();
		irrlichtVis->AddCamera(ChVector<>(0, -30, 0), ChVector<>(0, 0, 0));
		irrlichtVis->AddTypicalLights();

		// add play/pause button
		bool buttonPressed = false;
		MyActionReceiver receiver(irrlichtVis.get(), buttonPressed);
		irrlichtVis->AddUserEventReceiver(&receiver);

		// main simulation loop
		while (irrlichtVis->Run() && system.GetChTime() <= simulationDuration) {
			irrlichtVis->BeginScene();
			irrlichtVis->Render();
			irrlichtVis->EndScene();
			if (buttonPressed) {
				// step the simulation forwards
				system.DoStepDynamics(timestep);
				// append data to output file
				zPosition << std::left << std::setw(10) << std::setprecision(2) << std::fixed << system.GetChTime()
					<< std::right << std::setw(12) << std::setprecision(4) << std::fixed << sphereBody->GetPos().z()
					<< std::right << std::setw(18) << std::setprecision(4) << std::fixed << sphereBody->GetPos_dt().z()
					<< std::right << std::setw(18) << std::setprecision(2) << std::fixed << sphereBody->GetAppliedForce().z()
					<< std::endl;
				// force playback to be real-time
				// realtime_timer.Spin(timestep);
			}
		}
		zPosition.close();
		//auto end = std::chrono::high_resolution_clock::now();
		//unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		//std::cout << "Duration: " << duration/1000.0 << " seconds" << std::endl;
	}
	else{
		int frame = 0;
		while (system.GetChTime() <= simulationDuration) {
			system.DoStepDynamics(timestep);
			// append data to output file
			zPosition << std::left << std::setw(10) << std::setprecision(2) << std::fixed << system.GetChTime()
				<< std::right << std::setw(12) << std::setprecision(4) << std::fixed << sphereBody->GetPos().z()
				<< std::right << std::setw(18) << std::setprecision(4) << std::fixed << sphereBody->GetPos_dt().z()
				<< std::right << std::setw(18) << std::setprecision(2) << std::fixed << sphereBody->GetAppliedForce().z()
				<< std::endl;
			zPosition.close();
			frame++;
		}
	}
	return 0;
}