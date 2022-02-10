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
	ChIrrApp application(&system, L"Sphere Decay Test", core::dimension2d<u32>(800, 600), VerticalDir::Z);
	application.AddTypicalLogo();
	application.AddTypicalSky();
	application.AddTypicalLights();
	application.AddTypicalCamera(core::vector3df(-7.5, 30, 0), core::vector3df(-7.5, 0, 0)); // arguments are (location, orientation) as vectors

	// make first body with no force (other than gravity)
	auto body_1 = chrono_types::make_shared<ChBodyEasySphere>(5, 1); // arguments are (radius, density) TODO: what density to use?
	system.AddBody(body_1);
	body_1->SetPos(ChVector<>(-15, 0, -1));
	body_1->SetIdentifier(1);
	body_1->SetBodyFixed(false);
	body_1->SetCollide(false);
	body_1->SetMass(1);
	// Note ChBody::SetInertia can be used to set the entire Inertia matrix
	body_1->SetInertiaXX(ChVector<>(1, 1, 1));

	// Attach a visualization asset.
	auto sph_1 = chrono_types::make_shared<ChSphereShape>();
	body_1->AddAsset(sph_1);
	auto col_1 = chrono_types::make_shared<ChColorAsset>();
	col_1->SetColor(ChColor(0.6f, 0, 0));
	body_1->AddAsset(col_1);

	// set up body_2 (with forces in addition to gravity) with a mesh
	std::shared_ptr<ChBody> body_2 = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../test_for_chrono/oes_task10_sphere.obj").c_str(),                 // file name
		1000,                                                                                     // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false,                                                                                    // do not collide
		nullptr,                                                                                  // no need for contact material
		0                                                                                         // swept sphere radius
		);
	system.Add(body_2);
	body_2->SetPos(ChVector<>(0, 0, -1));
	body_2->SetMass(261.8e3);
	// attach color asset to body_2
	auto col_2 = chrono_types::make_shared<ChColorAsset>();
	col_2->SetColor(ChColor(0, 0, 0.6f));
	body_2->AddAsset(col_2);

	// testing adding external forces to the body_2
	BodyFileInfo sphere_file_info("../../test_for_chrono/sphere.h5", "body1");
	LinRestorForce lin_restor_force_2(sphere_file_info, body_2);
	ImpulseResponseForce irf(sphere_file_info, body_2);
	// declare some forces to be initialized in lin_restor_force_2 to be applied to to body_2 later
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
	body_2->AddForce(force);
	body_2->AddForce(torque);
	body_2->AddForce(force2);
	body_2->AddForce(torque2);

	// add buoyancy force from h5 file info
	auto fb = chrono_types::make_shared<BuoyancyForce>(sphere_file_info);
	body_2->AddForce(fb->getForce_ptr());

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

	// set up output file for body_2 position each step
	std::ofstream zpos("outfile/output.txt", std::ofstream::out);
	zpos.precision(5);
	zpos.width(10);
	//zpos.SetNumFormat("%10.5f"); // 10 characters displayed total, 5 digit precision (after decimal)

	// Simulation loop
	int frame = 0;
	//bool full_period = false;
	//ChVector<> initial_pos = body_2->GetPos();
	zpos << "#Time\t\tBody_2 Pos\n";
	while (application.GetDevice()->run()) {
		application.BeginScene();

		application.DrawAll();
		if (buttonPressed) {
			application.DoStep();
			zpos << system.GetChTime() << "\t" << body_2->GetPos().z() << /*"\t" << body_2->GetPos_dt().z() <<*/ "\n";
			//std::cout << body_2->GetChTime() << "\t" << body_2->GetPos()[2] << "\t" << body_2->GetPos_dt()[2] << "\t" << body_2->GetAppliedForce()[2] << "\n";
			frame++;
			//if (!full_period && (body_2->GetPos().Equals(initial_pos, 0.00001) ) && frame > 5 ) {
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
