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
	system.Set_G_acc(ChVector<>(0, 0, -9.81));
	double timestep = 0.06;
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
	std::vector<double> plate_heave_position;

	// set up body from a mesh
	if (!std::filesystem::exists("../../HydroChrono/demos/rm3/geometry/float.obj")) {
		std::cout << "File " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/demos/rm3/geometry/float.obj").c_str()) << " does not exist" << std::endl;
		return 0;
	}
	//std::cout << "Attempting to open mesh file: " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/meshFiles/float.obj").c_str()) << std::endl;
	std::shared_ptr<ChBody> float_body1 = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/demos/rm3/geometry/float.obj").c_str(),                 // file name
		1000,                                                                                     // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false,                                                                                    // do not collide
		nullptr,                                                                                  // no need for contact material
		0                                                                                         // swept sphere radius
		);

	// set up body from a mesh
	if (!std::filesystem::exists("../../HydroChrono/demos/rm3/geometry/plate.obj")) {
		std::cout << "File " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/demos/rm3/geometry/plate.obj").c_str()) << " does not exist" << std::endl;
		return 0;
	}
	//std::cout << "Attempting to open mesh file: " << std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/meshFiles/plate.obj").c_str()) << std::endl;
	std::shared_ptr<ChBody> plate_body2 = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/demos/rm3/geometry/plate.obj").c_str(),                 // file name
		1000,                                                                                     // density
		false,                                                                                                                                                                                // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false,                                                                                    // do not collide
		nullptr,                                                                                  // no need for contact material
		0                                                                                         // swept sphere radius
		);

	// define the float's initial conditions
	system.Add(float_body1);
	float_body1->SetNameString("body1"); 
	//float_body1->SetPos(ChVector<>(0, 0, 0));
	float_body1->SetMass(725.834);
	float_body1->SetCollide(false);

	// define the plate's initial conditions
	system.Add(plate_body2);
	plate_body2->SetNameString("body2");
	//plate_body2->SetPos(ChVector<>(0, 0, 0));
	plate_body2->SetMass(886.691);
	plate_body2->SetCollide(false);

	// TODO: add constraint so they only move up and down!

	// define wave parameters (not used in this demo)
	HydroInputs my_hydro_inputs;
	my_hydro_inputs.SetRegularWaveAmplitude(0.022);
	my_hydro_inputs.SetRegularWaveOmega(2.10);

	// attach hydrodynamic forces to body
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(float_body1);
	bodies.push_back(plate_body2);
	TestHydro blah(bodies, "../../HydroChrono/demos/rm3/hydroData/rm3.h5", my_hydro_inputs);

	// for profiling
	auto start = std::chrono::high_resolution_clock::now();

	if (visualizationOn) {
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
				float_heave_position.push_back(float_body1->GetPos().z());
				plate_heave_position.push_back(plate_body2->GetPos().z());
				// force playback to be real-time
				// realtime_timer.Spin(timestep);
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
			float_heave_position.push_back(float_body1->GetPos().z());
			plate_heave_position.push_back(plate_body2->GetPos().z());
			frame++;
		}
	}

	// for profiling
	auto end = std::chrono::high_resolution_clock::now();
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	if (profilingOn) {
		std::ofstream profilingFile;
		profilingFile.open("./results/rm3/decay/duration_ms.txt");
		profilingFile << duration << "\n";
		profilingFile.close();
	}

	if (saveDataOn) {
		std::ofstream outputFile;
		outputFile.open("./results/rm3/decay/rm3_decay.txt");
		if (!outputFile.is_open()) {
			std::cout << "Cannot open file " << std::filesystem::absolute("./results/rm3/decay/rm3_decay.txt") << std::endl;
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



	//std::shared_ptr<ChLoadContainer> my_loadcontainer;
	//std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
	//my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();
	//ChMatrixDynamic<> amMatrix;

	//my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(amMatrix, bodies);

	//system.Add(my_loadcontainer);
	//my_loadcontainer->Add(my_loadbodyinertia);


	//// update irrlicht app with body info
	//application.AssetBindAll();
	//application.AssetUpdateAll();

	//// some tools to handle the pause button
	//bool buttonPressed = false;
	//MyEventReceiver receiver(&application, buttonPressed);
	//application.SetUserEventReceiver(&receiver);

	// Info about which solver to use - may want to change this later

	//auto gmres_solver = chrono_types::make_shared<ChSolverGMRES>();  // change to mkl or minres?
	//gmres_solver->SetMaxIterations(300);
	//system.SetSolver(gmres_solver);
	//double timestep = 0.06; // also sets the timesteps in chrono system
	//application.SetTimestep(timestep);

	//MyEventReceiver receiver(application.get());
	// note how to add the custom event receiver to the default interface:
	//application->AddUserEventReceiver(&receiver);

	//// Solver settings
	//system.SetSolverType(ChSolver::Type::GMRES);
	//system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
	//system.SetStep(0.06);

	//// Simulation loop
	//double timestep = 0.03;
	////ChRealtimeStepTimer realtime_timer;
	//while (application->Run() && system.GetChTime() <= 40) {
	//	application->BeginScene();
	//	application->Render();
	//	tools::drawGrid(application.get(), 2, 2, 30, 30, ChCoordsys<>(ChVector<>(0, 0.01, 0), Q_from_AngX(CH_C_PI_2)),
	//		ChColor(0.3f, 0.3f, 0.3f), true);

	//	application->EndScene();
	//	system.DoStepDynamics(timestep);
	//	//realtime_timer.Spin(timestep);
	//}

	//// TODO compare solvers?
	////auto gmres_solver = chrono_types::make_shared<ChSolverGMRES>();  // change to mkl or minres?
	////gmres_solver->SetMaxIterations(300);
	////system.SetSolver(gmres_solver);
	////double timestep = 0.06; // also sets the timesteps in chrono system
	////application.SetTimestep(timestep);


	//// set up output file for body position each step
	//std::string of = "output.txt";                    /// < put name of your output file here
	//std::ofstream zpos(of, std::ofstream::out);
	//if (!zpos.is_open()) {
	//	std::cout << "Error opening file \"" + of + "\". Please make sure this file path exists then try again\n";
	//	return -1;
	//}	
	//std::cout << "Writing positions to file: " << std::filesystem::absolute(of) << std::endl;
	//zpos.precision(10);
	//zpos.width(12);
	//zpos << "#Time\tHeave " << float_body1->GetNameString() << "\t" << plate_body2->GetNameString() << std::endl;
	////zpos << "#Time\tBody vel" << std::endl;


	//// Simulation loop
	//int frame = 0;

	//while (application->Run() && system.GetChTime() < 3) {
	//	application->BeginScene();
	//	application->Render();
	//	tools::drawAllCOGs(application.get(), 15); // draws all cog axis lines, kinda neat
	//	//tools::drawGrid(application.GetVideoDriver(), 4, 4);
	//	/*if (buttonPressed)*/if(true) {
	//		zpos << system.GetChTime() << "\t" << float_body1->GetPos_dt().x() << "\t" << float_body1->GetPos_dt().y() << "\t" << float_body1->GetPos_dt().z();
	//		zpos << "\t" << float_body1->GetWvel_par().x() << "\t" << float_body1->GetWvel_par().y() << "\t" << float_body1->GetWvel_par().z() << std::endl;
	//		zpos << system.GetChTime() << "\t" << plate_body2->GetPos_dt().x() << "\t" << plate_body2->GetPos_dt().y() << "\t" << plate_body2->GetPos_dt().z();
	//		zpos << "\t" << plate_body2->GetWvel_par().x() << "\t" << plate_body2->GetWvel_par().y() << "\t" << plate_body2->GetWvel_par().z() << std::endl;
	//		system.DoStepDynamics(timestep);

	//		frame++;
	//	}
	//	application->EndScene();
	//}
	//zpos.close();
	//auto end = std::chrono::high_resolution_clock::now();
	//unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	//std::cout << "Duration: " << duration/1000.0 << " seconds" << std::endl;
	return 0;
}