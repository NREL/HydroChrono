// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban, Zuriah Quinton
// =============================================================================
//
// Recall that Irrlicht uses a left-hand frame, so everything is rendered with
// left and right flipped.
//
// =============================================================================

#include "H5_force_classes.h"
#include "chrono_irrlicht/ChIrrNodeAsset.h"

using namespace irr;
using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;

// =============================================================================
class MyEventReceiver : public IEventReceiver {
public:
	MyEventReceiver(ChIrrAppInterface* myapp, bool& buttonPressed)
		: pressed(buttonPressed) {
		// store pointer application
		application = myapp;

		// ..add a GUI button to control pause/play
		pauseButton = application->GetIGUIEnvironment()->addButton(rect<s32>(510, 20, 650, 35));
		buttonText = application->GetIGUIEnvironment()->addStaticText(L"Paused", rect<s32>(560, 20, 600, 35), false);
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
	ChIrrAppInterface* application;
	IGUIButton* pauseButton;
	IGUIStaticText* buttonText;

	bool& pressed;
};

int main(int argc, char* argv[]) {
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0, 0, -9.81));

	// Create the Irrlicht application for visualizing
	// names of some functions changed feb 7, 2022 (AddTypicalLogo -> AddLogo, etc)
	ChIrrApp application(&system, L"Sphere Decay Test", core::dimension2d<u32>(800, 600), VerticalDir::Z);
	application.AddLogo();
	application.AddSkyBox();
	application.AddTypicalLights();
	application.AddCamera(core::vector3df(0, 30, 0), core::vector3df(0, 0, 0)); // arguments are (location, orientation) as vectors

	// set up body from a mesh
	std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../test_for_chrono/oes_task10_sphere.obj").c_str(),                 // file name
		1000,                                                                                     // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false,                                                                                    // do not collide
		nullptr,                                                                                  // no need for contact material
		0                                                                                         // swept sphere radius
		);
	
	// old shpere stuff (not mesh)
	//std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasySphere>(5, 1);
	//auto sph = chrono_types::make_shared<ChSphereShape>();
	//body->AddAsset(sph);

	system.Add(body);
	body->SetPos(ChVector<>(0, 0, -1));
	//body->SetMass(261.8e3);
	//body->SetMass(261.8e3 + 130.8768);
	body->SetMass(261.8e3 + 130.8768e3); // added mass 130.8768 kg times rho=1000
	body->SetInertiaXX(ChVector<>(1, 1, 1));
	// attach color asset to body
	auto col_2 = chrono_types::make_shared<ChColorAsset>();
	col_2->SetColor(ChColor(0, 0, 0.6f));
	body->AddAsset(col_2);

	// testing adding external forces to the body
	BodyFileInfo sphere_file_info("../../test_for_chrono/sphere.h5", "body1");
	LinRestorForce lin_restor_force_2(sphere_file_info, body);
	ImpulseResponseForce irf(sphere_file_info, body);
	// declare some forces to be initialized in lin_restor_force_2 to be applied to to body later
	auto force = chrono_types::make_shared<ChForce>();
	auto torque = chrono_types::make_shared<ChForce>();
	auto force2 = chrono_types::make_shared<ChForce>();
	auto torque2 = chrono_types::make_shared<ChForce>();
	// set torque flag for torque
	torque->SetMode(ChForce::ForceType::TORQUE);
	torque2->SetMode(ChForce::ForceType::TORQUE);
	// initialize force and torque with member functions
	lin_restor_force_2.SetForce(force);
	lin_restor_force_2.SetTorque(torque);
	irf.SetForce(force2);
	irf.SetTorque(torque2);
	// apply force and torque to the body
	body->AddForce(force);
	body->AddForce(torque);
	body->AddForce(force2);
	body->AddForce(torque2);

	// add buoyancy force from h5 file info
	auto fb = chrono_types::make_shared<BuoyancyForce>(sphere_file_info);
	body->AddForce(fb->getForce_ptr());

	// update irrlicht app with body info
	//application.AssetBind(body_1);
	application.AssetBindAll();
	application.AssetUpdateAll();

	// some tools to handle the pause button
	bool buttonPressed = false;
	MyEventReceiver receiver(&application, buttonPressed);
	application.SetUserEventReceiver(&receiver);

	// Info about which solver to use - may want to change this later
	auto gmres_solver = chrono_types::make_shared<ChSolverMINRES>();  // change to mkl or minres?
	gmres_solver->SetMaxIterations(300);
	system.SetSolver(gmres_solver);
	double timestep = 0.015; // also sets the timesteps in system it seems
	application.SetTimestep(timestep);

	// set up output file for body position each step
	std::ofstream zpos("outfile/output.txt", std::ofstream::out);
	zpos.precision(10);
	zpos.width(12);
	//zpos.SetNumFormat("%10.5f"); // 10 characters displayed total, 5 digit precision (after decimal)

	// Simulation loop
	int frame = 0;
	//bool full_period = false;
	//ChVector<> initial_pos = body->GetPos();

	std::cout << "Body mass=" << body->GetMass() << std::endl;
	zpos << "#Time\tBody Pos\tBody vel (heave)\tforce (heave)\n";
	system.EnableSolverMatrixWrite(true, "blah");
	while (application.GetDevice()->run() && system.GetChTime() <= 20) {
		application.BeginScene();
		application.DrawAll();
		/*if (buttonPressed)*/if(true) {
			if (frame == 8) {
				ChSparseMatrix M;
				system.GetMassMatrix(&M);
				std::cout << "initial mass matrix\n" << M << std::endl;
				system.DumpSystemMatrices(true, true, true, true, "C:\\Users\\ZQUINTON\\code\\test_for_chrono_build\\Release\\outfile\\decay");
			}
			zpos << system.GetChTime() << "\t" << body->GetPos().z() << "\t" << body->GetPos_dt().z() << "\t" << body->GetAppliedForce().z() << "\n";
			application.DoStep();
			frame++;
			//if (!full_period && (body->GetPos().Equals(initial_pos, 0.00001) ) && frame > 5 ) {
			//	full_period = true;
			//	std::cout << "frame: " << frame << std::endl;
			//}
			//GetLog() << "\n" << fb->getForce_ptr()->GetForce() << "that's the force pointer value\n\n";
		}

		application.EndScene();
	}
	zpos.close();
	return 0;
}
