#include "../../src/hydro_forces.h"
//#include "./src/hydro_forces.h"
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

int main(int argc, char* argv[]) {
	//auto start = std::chrono::high_resolution_clock::now();
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	// system/solver settings
	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
	double timestep = 0.02;
	system.SetTimestepperType(ChTimestepper::Type::HHT);
	system.SetSolverType(ChSolver::Type::GMRES);
	system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
	system.SetStep(timestep);
	ChRealtimeStepTimer realtime_timer;
	double simulationDuration = 3.0;

	// some io/viz options
	bool visualizationOn = true;
	bool profilingOn = true;
	bool saveDataOn = true;
	std::vector<double> time_vector;
	std::vector<double> float_heave_position;
	std::vector<double> plate_heave_position;

	// set up body from a mesh
	if (!std::filesystem::exists("../../HydroChrono/demos/oswec/geometry/flap.obj")) {
		std::cout << "File " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/demos/oswec/geometry/flap.obj").c_str()) << " does not exist" << std::endl;
		return 0;
	}
	//std::cout << "Attempting to open mesh file: " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/meshFiles/float.obj").c_str()) << std::endl;
	std::shared_ptr<ChBody> flap_body = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/demos/oswec/geometry/flap.obj").c_str(),                 // file name
		0,                                                                                        // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false                                                                                    // collisions
		);

	// set up body from a mesh
	if (!std::filesystem::exists("../../HydroChrono/demos/oswec/geometry/base.obj")) {
		std::cout << "File " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/demos/oswec/geometry/base.obj").c_str()) << " does not exist" << std::endl;
		return 0;
	}

	//std::cout << "Attempting to open mesh file: " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/meshFiles/plate.obj").c_str()) << std::endl;
	std::shared_ptr<ChBody> base_body = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/demos/oswec/geometry/base.obj").c_str(),                 // file name
		0,                                                                                        // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false                                                                                    // collisions
		);

	// define the float's initial conditions
	system.Add(flap_body);
	flap_body->SetNameString("body1");
	// flap_body->SetPos(ChVector<>(0, 0, -3.9));
	flap_body->SetPos(ChVector<>(1.05925388, 0., -3.99267271));
	flap_body->SetRot(Q_from_AngAxis(CH_C_PI / 18, VECT_Y));
	flap_body->SetMass(127000);
	flap_body->SetInertiaXX(ChVector<>(1.85e6, 1.85e6, 1.85e6));

	// define the plate's initial conditions
	system.Add(base_body);
	base_body->SetNameString("body2");
	base_body->SetPos(ChVector<>(0, 0, -10.9));
	base_body->SetMass(999);
	base_body->SetInertiaXX(ChVector<>(1, 1, 1));
	base_body->SetBodyFixed(true);

	// define base-fore flap joint
	ChVector<> revolutePos(0.0, 0.0, -10.0);
	//ChQuaternion<> revoluteRot = Q_from_AngX(0.0); // CH_C_PI / 2.0);
	auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
	revolute->Initialize(base_body, flap_body, ChCoordsys<>(revolutePos));// , revoluteForeRot));
	system.AddLink(revolute);

	// add rotational spring-damper:

	// define wave parameters (not used in this demo)
	HydroInputs my_hydro_inputs;
	my_hydro_inputs.mode = noWaveCIC;// or 'regular' or 'regularCIC' or 'irregular';
	//my_hydro_inputs.regular_wave_amplitude = 0.022;
	//my_hydro_inputs.regular_wave_omega = 2.10;

	// attach hydrodynamic forces to body
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(flap_body);
	bodies.push_back(base_body);
	TestHydro blah(bodies, "../../HydroChrono/demos/oswec/hydroData/oswec.h5", my_hydro_inputs);

	//// Debug printing added mass matrix and system mass matrix
	//ChSparseMatrix M;
	//system.GetMassMatrix(&M);
	//std::cout << M << std::endl;

	// for profiling
	auto start = std::chrono::high_resolution_clock::now();

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
		irrlichtVis->AddCamera(ChVector<>(0, -50, -10), ChVector<>(0, 0, -10)); // camera position and where it points
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
				//system.GetMassMatrix(&M);
				//std::cout << M << std::endl;
				// step the simulation forwards
				system.DoStepDynamics(timestep);
				// append data to std vector
				time_vector.push_back(system.GetChTime());
				float_heave_position.push_back(flap_body->GetPos().z());
				plate_heave_position.push_back(base_body->GetPos().z());
				// force playback to be real-time
				realtime_timer.Spin(timestep);
			}
		}
	}
	else {
		int frame = 0;
		while (system.GetChTime() <= simulationDuration) {
			// append data to std vector
			time_vector.push_back(system.GetChTime());
			float_heave_position.push_back(flap_body->GetPos().z());
			plate_heave_position.push_back(base_body->GetPos().z());

			// step the simulation forwards
			system.DoStepDynamics(timestep);

			frame++;
		}
	}

	// for profiling
	auto end = std::chrono::high_resolution_clock::now();
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	if (profilingOn) {
		std::ofstream profilingFile;
		profilingFile.open("./results/oswec/decay/duration_ms.txt");
		if (!profilingFile.is_open()) {
			if (!std::filesystem::exists("./results/oswec/decay")) {
				std::cout << "Path " << std::filesystem::absolute("./results/oswec/decay") << " does not exist, creating it now..." << std::endl;
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
				std::cout << "Path " << std::filesystem::absolute("./results/oswec/decay") << " does not exist, creating it now..." << std::endl;
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
		outputFile << std::left << std::setw(10) << "Time (s)"
			<< std::right << std::setw(16) << "Float Heave (m)"
			<< std::right << std::setw(16) << "Plate Heave (m)"
			<< std::endl;
		for (int i = 0; i < time_vector.size(); ++i)
			outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
			<< std::right << std::setw(16) << std::setprecision(4) << std::fixed << float_heave_position[i]
			<< std::right << std::setw(16) << std::setprecision(4) << std::fixed << plate_heave_position[i]
			<< std::endl;
		outputFile.close();
	}
	return 0;
}