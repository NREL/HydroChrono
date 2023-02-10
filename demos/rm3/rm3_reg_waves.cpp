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
	double timestep = 0.01;
	system.SetTimestepperType(ChTimestepper::Type::HHT);
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
	std::vector<double> float_heave_position;
	std::vector<double> float_drift_position;
	std::vector<double> plate_heave_position;

	// set up body from a mesh
	if (!std::filesystem::exists("../../HydroChrono/demos/rm3/geometry/float_cog.obj")) {
		std::cout << "File " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/demos/rm3/geometry/float.obj").c_str()) << " does not exist" << std::endl;
		return 0;
	}
	//std::cout << "Attempting to open mesh file: " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/meshFiles/float.obj").c_str()) << std::endl;
	std::shared_ptr<ChBody> float_body1 = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/demos/rm3/geometry/float_cog.obj").c_str(),                 // file name
		0,                                                                                        // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false                                                                                     // collisions
		);

	// set up body from a mesh
	if (!std::filesystem::exists("../../HydroChrono/demos/rm3/geometry/plate_cog.obj")) {
		std::cout << "File " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/demos/rm3/geometry/plate.obj").c_str()) << " does not exist" << std::endl;
		return 0;
	}

	//std::cout << "Attempting to open mesh file: " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/meshFiles/plate.obj").c_str()) << std::endl;
	std::shared_ptr<ChBody> plate_body2 = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/demos/rm3/geometry/plate_cog.obj").c_str(),                 // file name
		0,                                                                                        // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false                                                                                     // collisions
		);

	// define the float's initial conditions
	system.Add(float_body1);
	float_body1->SetNameString("body1"); 
	float_body1->SetPos(ChVector<>(0, 0, -0.72));
	float_body1->SetMass(725834);
	float_body1->SetInertiaXX(ChVector<>(20907301.0, 21306090.66, 37085481.11));
	//float_body1->SetCollide(false);

	// define the plate's initial conditions
	system.Add(plate_body2);
	plate_body2->SetNameString("body2");
	plate_body2->SetPos(ChVector<>(0, 0, (-21.29)));
	plate_body2->SetMass(886691);
	plate_body2->SetInertiaXX(ChVector<>(94419614.57, 94407091.24, 28542224.82));
	//plate_body2->SetCollide(false);

	// add prismatic joint between the two bodies
	auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
	prismatic->Initialize(float_body1, plate_body2, false, ChCoordsys<>(ChVector<>(0, 0, -0.72)), ChCoordsys<>(ChVector<>(0, 0, -21.29)));
	system.AddLink(prismatic);

	auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
	prismatic_pto->Initialize(float_body1, plate_body2, false, ChVector<>(0, 0, -0.72), ChVector<>(0, 0, -21.29));
	prismatic_pto->SetDampingCoefficient(0.0);
	system.AddLink(prismatic_pto);

	// define wave parameters 
	HydroInputs my_hydro_inputs;
	my_hydro_inputs.mode = regular;
	my_hydro_inputs.regular_wave_amplitude = 1.0;
	my_hydro_inputs.regular_wave_omega = 2.10;

	// attach hydrodynamic forces to body
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(float_body1);
	bodies.push_back(plate_body2);
	TestHydro blah(bodies, "../../HydroChrono/demos/rm3/hydroData/rm3.h5", my_hydro_inputs);

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
		irrlichtVis->SetWindowTitle("RM3 - Regular Wave Test");
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
				// step the simulation forwards
				system.DoStepDynamics(timestep);
				// append data to std vector
				time_vector.push_back(system.GetChTime());
				float_heave_position.push_back(float_body1->GetPos().z());
				float_drift_position.push_back(float_body1->GetPos().x());
				plate_heave_position.push_back(plate_body2->GetPos().z());
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
			float_heave_position.push_back(float_body1->GetPos().z());
			float_drift_position.push_back(float_body1->GetPos().x());
			plate_heave_position.push_back(plate_body2->GetPos().z());

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
		profilingFile.open("./results/rm3/reg_waves/duration_ms.txt");
		if (!profilingFile.is_open()) {
			if (!std::filesystem::exists("./results/rm3/reg_waves")) {
				std::cout << "Path " << std::filesystem::absolute("./results/rm3/reg_waves") << " does not exist, creating it now..." << std::endl;
				std::filesystem::create_directory("./results");
				std::filesystem::create_directory("./results/rm3");
				std::filesystem::create_directory("./results/rm3/reg_waves");
				profilingFile.open("./results/rm3/reg_waves/duration_ms.txt");
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
		outputFile.open("./results/rm3/reg_waves/rm3_reg_waves.txt");
		if (!outputFile.is_open()) {
			if (!std::filesystem::exists("./results/rm3/reg_waves")) {
				std::cout << "Path " << std::filesystem::absolute("./results/rm3/reg_waves") << " does not exist, creating it now..." << std::endl;
				std::filesystem::create_directory("./results");
				std::filesystem::create_directory("./results/rm3");
				std::filesystem::create_directory("./results/rm3/reg_waves");
				outputFile.open("./results/rm3/reg_waves/rm3_decay.txt");
				if (!outputFile.is_open()) {
					std::cout << "Still cannot open file, ending program" << std::endl;
					return 0;
				}
			}
		}
		outputFile << std::left << std::setw(10) << "Time (s)"
			<< std::right << std::setw(16) << "Float Heave (m)"
			<< std::right << std::setw(16) << "Plate Heave (m)"
			<< std::right << std::setw(16) << "Float Drift (x) (m)"
			<< std::endl;
		for (int i = 0; i < time_vector.size(); ++i)
			outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
			<< std::right << std::setw(16) << std::setprecision(4) << std::fixed << float_heave_position[i]
			<< std::right << std::setw(16) << std::setprecision(4) << std::fixed << plate_heave_position[i]
			<< std::right << std::setw(16) << std::setprecision(4) << std::fixed << float_drift_position[i]
			<< std::endl;
		outputFile.close();
	}
	return 0;
}