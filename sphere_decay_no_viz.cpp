#include "hydro_forces.h"
#include "chrono_irrlicht/ChIrrNodeAsset.h"
#include <chrono>

//using namespace irr;
//using namespace irr::core;
//using namespace irr::scene;
//using namespace irr::video;
//using namespace irr::io;
//using namespace irr::gui;

// =============================================================================
//class MyEventReceiver : public IEventReceiver {
//public:
//	MyEventReceiver(ChIrrAppInterface* myapp, bool& buttonPressed)
//		: pressed(buttonPressed) {
//		// store pointer application
//		application = myapp;
//
//		// ..add a GUI button to control pause/play
//		pauseButton = application->GetIGUIEnvironment()->addButton(rect<s32>(510, 20, 650, 35));
//		buttonText = application->GetIGUIEnvironment()->addStaticText(L"Paused", rect<s32>(560, 20, 600, 35), false);
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
//	ChIrrAppInterface* application;
//	IGUIButton* pauseButton;
//	IGUIStaticText* buttonText;
//
//	bool& pressed;
//};

int main(int argc, char* argv[]) {
	auto start = std::chrono::high_resolution_clock::now();
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0, 0, -9.81));

	// Create the Irrlicht application for visualizing
	//ChIrrApp application(&system, L"Sphere Decay Test", core::dimension2d<u32>(800, 600), VerticalDir::Z);
	//application.AddLogo();
	//application.AddSkyBox();
	//application.AddTypicalLights();
	//application.AddCamera(core::vector3df(0, 30, 0), core::vector3df(0, 0, 0)); // arguments are (location, orientation) as vectors

	// set up body from a mesh
	//std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasyMesh>(                   //
	//	GetChronoDataFile("../../test_for_chrono/oes_task10_sphere.obj").c_str(),                 // file name
	//	1000,                                                                                     // density
	//	false,                                                                                    // do not evaluate mass automatically
	//	true,                                                                                     // create visualization asset
	//	false,                                                                                    // do not collide
	//	nullptr,                                                                                  // no need for contact material
	//	0                                                                                         // swept sphere radius
	//	);
	
	// old sphere stuff (for when you're not using mesh above)
	std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasySphere>(5, 1);
	auto sph = chrono_types::make_shared<ChSphereShape>();
	body->AddAsset(sph);

	// set up body initial conditions
	system.Add(body);
	body->SetPos(ChVector<>(0, 0, -2));
	body->SetMass(261.8e3);
	// attach color asset to body
	auto col_2 = chrono_types::make_shared<ChColorAsset>();
	col_2->SetColor(ChColor(0, 0, 0.6f));
	body->AddAsset(col_2);

	HydroInputs myHydroInputs;
	myHydroInputs.regularWaveAmplitude = 0.022;
	myHydroInputs.regularWaveOmega = 2.10;
	LoadAllHydroForces blah(body, "../../HydroChrono/sphere.h5", myHydroInputs);

	//// testing adding hydro forces to the body-----------------------------------------------------------------
	//BodyFileInfo sphere_file_info("../../test_for_chrono/sphere.h5", "body1");     /// < object to read h5 file info
	//LinRestorForce lin_restor_force_2(sphere_file_info, body);                     /// < object for linear restoring force
	//ImpulseResponseForce irf(sphere_file_info, body);                              /// < object for impulse restoring force
	//
	//// declare some forces to be initialized in lin_restor_force_2 to be applied to to body later
	//auto force = chrono_types::make_shared<ChForce>();
	//auto torque = chrono_types::make_shared<ChForce>();
	//auto force2 = chrono_types::make_shared<ChForce>();
	//auto torque2 = chrono_types::make_shared<ChForce>();
	//// set torque flag for torque
	//torque->SetMode(ChForce::ForceType::TORQUE);
	//torque2->SetMode(ChForce::ForceType::TORQUE);
	//// initialize force and torque with member functions
	//lin_restor_force_2.SetForce(force);
	//lin_restor_force_2.SetTorque(torque);
	//irf.SetForce(force2);
	//irf.SetTorque(torque2);
	//// apply force and torque to the body
	//body->AddForce(force);
	//body->AddForce(torque);
	//body->AddForce(force2);
	//body->AddForce(torque2);

	//// add buoyancy force from h5 file info
	//auto fb = chrono_types::make_shared<BuoyancyForce>(sphere_file_info);
	//body->AddForce(fb->getForce_ptr());

	//// set up added mass as load
	//auto my_loadcontainer = chrono_types::make_shared< ChLoadContainer>();
	//system.Add(my_loadcontainer);
	//auto my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(body, sphere_file_info);
	//my_loadcontainer->Add(my_loadbodyinertia);
	// end hydro forces ----------------------------------------------------------------------------------

	// update irrlicht app with body info
	//application.AssetBindAll();
	//application.AssetUpdateAll();

	// some tools to handle the pause button
	//bool buttonPressed = false;
	//MyEventReceiver receiver(&application, buttonPressed);
	//application.SetUserEventReceiver(&receiver);

	// Info about which solver to use - may want to change this later
	auto gmres_solver = chrono_types::make_shared<ChSolverGMRES>();  // change to mkl or minres?
	gmres_solver->SetMaxIterations(300);
	system.SetSolver(gmres_solver);
	double timestep = 0.015; // also sets the timesteps in chrono system
	//system.SetTimestep(timestep);

	// set up output file for body position each step
	std::string of = "output.txt";                    /// < put name of your output file here
	std::ofstream zpos(of, std::ofstream::out);
	if (!zpos.is_open()) { 
		std::cout << "Error opening file \"" + of + "\". Please make sure this file path exists then try again\n";
		return -1;
	}
	zpos.precision(10);
	zpos.width(12);
	zpos << "#Time\tBody Pos\tBody vel (heave)\tforce (heave)\n";

	// Simulation loop
	int frame = 0;
	while (/*application.GetDevice()->run() && */system.GetChTime() <= 400) {
		//application.BeginScene();
		//application.DrawAll();
		/*if (buttonPressed)*/if(true) {
			zpos << system.GetChTime() << "\t" << body->GetPos().x() << "\t" << body->GetPos().z() << "\t" << body->GetPos_dt().z() << "\t" << body->GetAppliedForce().z() << "\n";
			system.DoStepDynamics(timestep);
			frame++;
		}
		/*application.EndScene();*/
	}
	zpos.close();
	auto end = std::chrono::high_resolution_clock::now();
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Duration: " << duration/1000.0 << " seconds" << std::endl;
	return 0;
}