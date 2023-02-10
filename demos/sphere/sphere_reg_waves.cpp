#include "./src/hydro_forces.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono_irrlicht/ChIrrMeshTools.h"
#include "chrono/core/ChRealtimeStep.h"
#include <filesystem>
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

	// Setup Ground
	auto ground = chrono_types::make_shared<ChBody>();
	system.AddBody(ground);
	ground->SetPos(ChVector<>(0, 0, -5));
	ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(false);

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
	sphereBody->SetPos(ChVector<>(0, 0, -2));
	sphereBody->SetMass(261.8e3);

	int reg_wave_num = 10;
	double task10_wave_amps[] = { 0.044, 0.078, 0.095, 0.123, 0.177, 0.24, 0.314, 0.397, 0.491, 0.594 };
	double task10_wave_omegas[] = { 2.094395102, 1.570796327, 1.427996661, 1.256637061, 1.047197551, 0.897597901, 0.785398163, 0.698131701, 0.628318531, 0.571198664 };
	double task10dampings[] = { 398736.034, 118149.758, 90080.857, 161048.558, 322292.419, 479668.979, 633979.761, 784083.286, 932117.647, 1077123.445 };
	//int waveNum = 0;

	// add prismatic joint between sphere and ground (limit to heave motion only)
	auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
	prismatic->Initialize(sphereBody, ground, false, ChCoordsys<>(ChVector<>(0, 0, -2)), ChCoordsys<>(ChVector<>(0, 0, -5)));
	system.AddLink(prismatic);

	// Create the spring between body_1 and ground. The spring end points are
	// specified in the body relative frames.
	double rest_length = 3.0;
	double spring_coef = 0.0;
	double damping_coef = task10dampings[reg_wave_num - 1];
	auto spring_1 = chrono_types::make_shared<ChLinkTSDA>();
	spring_1->Initialize(sphereBody, ground, false, ChVector<>(0, 0, -2), ChVector<>(0, 0, -5)); // false means positions are in global frame
	//spring_1->SetRestLength(rest_length); // if not set, the rest length is calculated from initial position
	spring_1->SetSpringCoefficient(spring_coef);
	spring_1->SetDampingCoefficient(damping_coef);
	system.AddLink(spring_1);

	HydroInputs my_hydro_inputs;
	//my_hydro_inputs.SetRegularWaveAmplitude(task10_wave_amps[reg_wave_num - 1]); //0.594 (for wave 10)
	//my_hydro_inputs.SetRegularWaveOmega(task10_wave_omegas[reg_wave_num - 1]); //0.571198664 (for wave 10)
	my_hydro_inputs.mode = regular; // uses regular wave mode
	my_hydro_inputs.regular_wave_amplitude = task10_wave_amps[reg_wave_num - 1]; //0.095;
	my_hydro_inputs.regular_wave_omega = task10_wave_omegas[reg_wave_num - 1];//1.427996661;

	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(sphereBody);
	TestHydro blah(bodies, "../../HydroChrono/demos/sphere/hydroData/sphere.h5", my_hydro_inputs);

	// for profiling
	auto start = std::chrono::high_resolution_clock::now();
	
	if (visualizationOn) {
		// create the irrlicht application for visualizing
		auto irrlichtVis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
		irrlichtVis->AttachSystem(&system);
		irrlichtVis->SetWindowSize(1280, 720);
		irrlichtVis->SetWindowTitle("Sphere - Regular Waves Test");
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

		// Main simulation loop
		int frame = 0;
		while (irrlichtVis->Run() && system.GetChTime() <= simulationDuration) {
			irrlichtVis->BeginScene();
			irrlichtVis->Render();
			irrlichtVis->EndScene();
			if (buttonPressed) {
				// step simulation forward
				system.DoStepDynamics(timestep);
				// append data to std vector
				time_vector.push_back(system.GetChTime());
				heave_position.push_back(sphereBody->GetPos().z());				
				frame++;
			}
		}
	}
	else {
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
		std::string out_dir = "results/regular_waves/";
		std::string out_file = "duration_ms.txt";
		std::ofstream profilingFile(out_dir + out_file);
		if (!profilingFile.is_open()) {
			if (!std::filesystem::exists(out_dir)) {
				std::cout << "Path " << std::filesystem::absolute(out_dir) << " does not exist, creating it now..." << std::endl;
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
		std::string out_dir = "results/regular_waves/";
		std::string out_file = "regwave_" + std::to_string(reg_wave_num) + ".txt";
		std::ofstream outputFile(out_dir + out_file);
		if (!outputFile.is_open()) {
			if (!std::filesystem::exists(out_dir)) {
				std::cout << "Path " << std::filesystem::absolute(out_dir) << " does not exist, creating it now..." << std::endl;
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
		outputFile << "Wave #: \t" << reg_wave_num << "\n";
		outputFile << "Wave amplitude (m): \t" << my_hydro_inputs.regular_wave_amplitude << "\n";
		outputFile << "Wave omega (rad/s): \t" << my_hydro_inputs.regular_wave_omega << "\n";
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
