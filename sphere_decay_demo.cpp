#include "hydro_forces.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono_irrlicht/ChIrrMeshTools.h"

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::irrlicht;

using namespace irr;
using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;

// =============================================================================
//class MyEventReceiver : public IEventReceiver {
//public:
//	MyEventReceiver(ChVisualSystemIrrlicht* myapp, bool& buttonPressed)
//		: pressed(buttonPressed) {
//		// store pointer application
//		application = myapp;
//
//		// ..add a GUI button to control pause/play
//		//pauseButton = application->GetIGUIEnvironment()->addButton(rect<s32>(510, 20, 650, 35));
//		//buttonText = application->GetIGUIEnvironment()->addStaticText(L"Paused", rect<s32>(560, 20, 600, 35), false);
//	}
//
//	bool OnEvent(const SEvent& event) {
//		// check if user clicked button
//		if (event.EventType == EET_GUI_EVENT) {
//			switch (event.GUIEvent.EventType) {
//			case EGET_BUTTON_CLICKED:
//				pressed = !pressed;
//				if (pressed) {
//					buttonText->setText(L"Playing");
//				}
//				else {
//					buttonText->setText(L"Paused");
//				}
//				return pressed;
//				break;
//			default:
//				break;
//			}
//		}
//		return false;
//	}
//
//private:
//	ChVisualSystemIrrlicht* application;
//	IGUIButton* pauseButton;
//	IGUIStaticText* buttonText;
//
//	bool& pressed;
//};

int main(int argc, char* argv[]) {
	//auto start = std::chrono::high_resolution_clock::now();
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";
	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0, 0, -9.81));

	//// Create the Irrlicht application for visualizing
	//ChIrrApp application(&system, L"Sphere Decay Test", core::dimension2d<u32>(800, 600), VerticalDir::Z);
	//application.AddLogo();
	//application.AddSkyBox();
	//application.AddTypicalLights();
	//application.AddCamera(core::vector3df(0, 30, 0), core::vector3df(0, 0, 0)); // arguments are (location, orientation) as vectors

	// Create the Irrlicht application for visualizing
	auto application = chrono_types::make_shared<ChVisualSystemIrrlicht>();

	application->AttachSystem(&system);
	application->SetWindowSize(800, 600);
	application->SetWindowTitle("Sphere Decay Test");
	application->SetCameraVertical(CameraVerticalDir::Y);
	application->Initialize();
	application->AddLogo();
	application->AddSkyBox();
	application->AddCamera(ChVector<>(0, 0, -6), ChVector<>(-2, 2, 0));
	application->AddTypicalLights();


	// set up body from a mesh
	std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/meshFiles/oes_task10_sphere.obj").c_str(),                 // file name
		1000,                                                                                     // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false,                                                                                    // do not collide
		nullptr,                                                                                  // no need for contact material
		0                                                                                         // swept sphere radius
		);

	// old sphere stuff (for when you're not using mesh above)
	//std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasySphere>(5, 1);
	//auto sph = chrono_types::make_shared<ChSphereShape>();
	//body->AddAsset(sph);

	// set up body initial conditions
	system.Add(body);
	body->SetNameString("body1"); // must set body name!
	body->SetPos(ChVector<>(0, 0, -1));
	body->SetMass(261.8e3);
	// attach color asset to body
	//auto col_2 = chrono_types::make_shared<ChColorAsset>();
	//col_2->SetColor(ChColor(0, 0, 0.6f));
	//body->AddAsset(col_2);

	HydroInputs my_hydro_inputs;
	my_hydro_inputs.SetRegularWaveAmplitude(0.022);
	my_hydro_inputs.SetRegularWaveOmega(2.10);
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(body);
	TestHydro blah(bodies, "../../HydroChrono/sphere.h5", my_hydro_inputs);
	// update irrlicht app with body info
	//application.AssetBindAll();
	//application.AssetUpdateAll();

	//MyEventReceiver receiver(application.get());
	// note how to add the custom event receiver to the default interface:
	//application->AddUserEventReceiver(&receiver);

	// Solver settings
	double timestep = 0.06;
	system.SetSolverType(ChSolver::Type::GMRES);
	system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
	system.SetStep(timestep);

	//// some tools to handle the pause button
	//bool buttonPressed = false;
	//MyEventReceiver receiver(&application, buttonPressed);
	//application.SetUserEventReceiver(&receiver);

	//// Info about which solver to use - may want to change this later
	//auto gmres_solver = chrono_types::make_shared<ChSolverGMRES>();  // change to mkl or minres?
	//gmres_solver->SetMaxIterations(300);
	//system.SetSolver(gmres_solver);
	//double timestep = 0.015; // also sets the timesteps in chrono system
	//application.SetTimestep(timestep);

	// set up output file for body position each step
	std::string of = "sphere_decay.txt";                    /// < put name of your output file here
	std::ofstream zpos(of, std::ofstream::out);
	if (!zpos.is_open()) {
		std::cout << "Error opening file \"" + of + "\". Please make sure this file path exists then try again\n";
		return -1;
	}
	zpos.precision(10);
	zpos.width(12);
	zpos << "#Time\tBody Pos\tBody vel (heave)\tforce (heave)" << std::endl;

	// Simulation loop
	int frame = 0;
	while (application->Run() && system.GetChTime() <= 40) {
		application->BeginScene();
		application->Render();
		/*if (buttonPressed)*/if(true) {
			zpos << system.GetChTime() << "\t" << body->GetPos().z() << "\t" << body->GetPos_dt().z() << "\t" << body->GetAppliedForce().z() << std::endl;
			//application.DoStep();
			system.DoStepDynamics(timestep);
			frame++;
		}
		application->EndScene();
	}
	zpos.close();
	//auto end = std::chrono::high_resolution_clock::now();
	//unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	//std::cout << "Duration: " << duration/1000.0 << " seconds" << std::endl;
	return 0;
}
