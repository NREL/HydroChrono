//#include "../../src/hydro_forces.h"
#include "./src/hydro_forces.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono_irrlicht/ChIrrMeshTools.h"
#include "chrono/core/ChRealtimeStep.h"
#include <iomanip> // std::setprecision
#include <chrono> // std::chrono::high_resolution_clock::now
#include <vector> // std::vector<double>

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
	GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";

	// system/solver settings
	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
	double timestep = 0.015;
	system.SetSolverType(ChSolver::Type::GMRES);
	system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
	system.SetStep(timestep);
	ChRealtimeStepTimer realtime_timer;
	double simulationDuration = 40.0;

	// some io/viz options
	bool visualizationOn = true;
	bool profilingOn = true;
	bool saveDataOn = true;
	std::vector<double> time_vector;
	std::vector<double> heave_position;

	// set up body from a mesh
	std::cout << "Attempting to open mesh file: " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/demos/sphere/geometry/oes_task10_sphere.obj").c_str()) << std::endl;
	std::shared_ptr<ChBody> sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(       //
		GetChronoDataFile("../../HydroChrono/demos/sphere/geometry/oes_task10_sphere.obj").c_str(),    // file name
		1000,                                                                             // density
		false,                                                                            // do not evaluate mass automatically
		true,                                                                             // create visualization asset
		false                                                                             // do not collide
		);

	// define the body's initial conditions
	system.Add(sphereBody);
	sphereBody->SetNameString("body1"); // must set body name correctly! (must match .h5 file)
	sphereBody->SetPos(ChVector<>(0, 0, -1));
	sphereBody->SetMass(261.8e3);

	// define wave parameters (not used in this demo)
	// Todo define a way to use TestHydro without hydro_inputs/waves
	HydroInputs my_hydro_inputs;
	my_hydro_inputs.mode = noWaveCIC;
	//my_hydro_inputs.regular_wave_amplitude = 0.022;
	//my_hydro_inputs.regular_wave_omega = 2.10;

	// attach hydrodynamic forces to body
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(sphereBody);
	TestHydro blah(bodies, "../../HydroChrono/demos/sphere/hydroData/sphere.h5", my_hydro_inputs);

	// for profiling
	auto start = std::chrono::high_resolution_clock::now(); 

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
				// append data to std vector
				time_vector.push_back(system.GetChTime());
				heave_position.push_back(sphereBody->GetPos().z());
				// force playback to be real-time
				// realtime_timer.Spin(timestep);
			}
		}
	}
	else{
		int frame = 0;
		while (system.GetChTime() <= simulationDuration) {
			// step the simulation forwards
			system.DoStepDynamics(timestep);
			// append data to std vector
			time_vector.push_back(system.GetChTime());
			heave_position.push_back(sphereBody->GetPos().z());
			frame++;
		}
	}

	// for profiling
	auto end = std::chrono::high_resolution_clock::now();
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	if (profilingOn) {
		std::ofstream profilingFile;
		profilingFile.open("./results/decay/duration_ms.txt");
		if (!profilingFile.is_open()) {
			if (!std::filesystem::exists("./results/decay")) {
				std::cout << "Path " << std::filesystem::absolute("./results/decay") << " does not exist, creating it now..." << std::endl;
				std::filesystem::create_directory("./results");
				std::filesystem::create_directory("./results/decay");
				profilingFile.open("./results/decay/duration_ms.txt");
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
		outputFile.open("./results/decay/sphere_decay.txt");
		if (!outputFile.is_open()) {
			if (!std::filesystem::exists("./results/decay")) {
				std::cout << "Path " << std::filesystem::absolute("./results/decay") << " does not exist, creating it now..." << std::endl;
				std::filesystem::create_directory("./results");
				std::filesystem::create_directory("./results/decay");
				outputFile.open("./results/decay/sphere_decay.txt");
				if (!outputFile.is_open()) {
					std::cout << "Still cannot open file, ending program" << std::endl;
					return 0;
				}
			}
		}
		outputFile << std::left << std::setw(10) << "Time (s)"
		<< std::right << std::setw(12) << "Heave (m)"
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