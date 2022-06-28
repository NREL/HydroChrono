#include "hydro_forces.h"
#include "chrono_irrlicht/ChIrrNodeAsset.h"
#include <chrono>

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
	auto start = std::chrono::high_resolution_clock::now();
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0, 0, -9.81));
	// Create the Irrlicht application for visualizing
	ChIrrApp application(&system, L"MultiBody Demo", core::dimension2d<u32>(800, 600), VerticalDir::Z);
	application.AddLogo();
	application.AddSkyBox();
	application.AddTypicalLights();
	application.AddCamera(core::vector3df(0, 70, -10), core::vector3df(0, 0, -10)); // arguments are (location, orientation) as vectors

	// uncomment out lines for xy plane on z=0 to be shown
	auto ground = chrono_types::make_shared<ChBodyEasyBox>(50, 50, 0.05, 10, true, false);
	system.AddBody(ground);
	ground->SetBodyFixed(true);
	ground->SetCollide(false);
	ground->SetPos(chrono::ChVector(0, 0, 0));


	// set up body from a mesh
	std::shared_ptr<ChBody> float_body1 = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/meshFiles/float.obj").c_str(),                 // file name
		1000,                                                                                     // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false,                                                                                    // do not collide
		nullptr,                                                                                  // no need for contact material
		0                                                                                         // swept sphere radius
		);

	// set up body from a mesh
	std::shared_ptr<ChBody> plate_body2 = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../HydroChrono/meshFiles/plate.obj").c_str(),                 // file name
		1000,                                                                                     // density
		false,                                                                                                                                                                                // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false,                                                                                    // do not collide
		nullptr,                                                                                  // no need for contact material
		0                                                                                         // swept sphere radius
		);

	// set up body initial conditions
	system.Add(float_body1);
	float_body1->SetNameString("body1"); // TODO do i want this?
	//float_body1->SetPos(ChVector<>(0, 0, 0));
	float_body1->SetMass(886.691 * 1000);
	float_body1->SetCollide(true);
	// attach color asset to body
	auto col_1 = chrono_types::make_shared<ChColorAsset>();
	col_1->SetColor(ChColor(0, 0, 0.6f));
	float_body1->AddAsset(col_1);

	// set up body initial conditions
	system.Add(plate_body2);
	plate_body2->SetNameString("body2");
	//plate_body2->SetPos(ChVector<>(0, 0, 0));
	plate_body2->SetMass(886.691 * 1000);
	plate_body2->SetCollide(true);
	// attach color asset to body
	auto col_2 = chrono_types::make_shared<ChColorAsset>();
	col_2->SetColor(ChColor(0, 0.7f, 0.8f));
	plate_body2->AddAsset(col_2);

	HydroInputs my_hydro_inputs;
	my_hydro_inputs.SetRegularWaveAmplitude(0.022);
	my_hydro_inputs.SetRegularWaveOmega(2.10);

	std::vector<std::shared_ptr<ChBody>> bodies = { float_body1, plate_body2 };
	//std::shared_ptr<ChLoadContainer> my_loadcontainer;
	//std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
	//my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();
	//ChMatrixDynamic<> amMatrix;

	//my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(amMatrix, bodies);

	//system.Add(my_loadcontainer);
	//my_loadcontainer->Add(my_loadbodyinertia);
	// TODO figure out better way for multibody
	// TODO if body's name is correct, then we dont need the name field here...not sure if thats how we want to do this
	LoadAllHydroForces hydroForcesTorus(float_body1, "../../HydroChrono/rm3.h5", float_body1->GetNameString(), my_hydro_inputs);
	LoadAllHydroForces hydroForcesCylinder(plate_body2, "../../HydroChrono/rm3.h5", plate_body2->GetNameString(), my_hydro_inputs);


	// update irrlicht app with body info
	application.AssetBindAll();
	application.AssetUpdateAll();

	// some tools to handle the pause button
	bool buttonPressed = false;
	MyEventReceiver receiver(&application, buttonPressed);
	application.SetUserEventReceiver(&receiver);

	// Info about which solver to use - may want to change this later
	auto gmres_solver = chrono_types::make_shared<ChSolverGMRES>();  // change to mkl or minres?
	gmres_solver->SetMaxIterations(300);
	system.SetSolver(gmres_solver);
	double timestep = 0.06; // also sets the timesteps in chrono system
	application.SetTimestep(timestep);

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
	while (application.GetDevice()->run() && system.GetChTime() < 1) {
		application.BeginScene();
		application.DrawAll();
		tools::drawAllCOGs(system, application.GetVideoDriver(), 15); // draws all cog axis lines, kinda neat
		//tools::drawGrid(application.GetVideoDriver(), 4, 4);
		/*if (buttonPressed)*/if(true) {
			zpos << system.GetChTime() << "\t" << float_body1->GetPos().z() << "\t" << float_body1->GetPos_dt().z() << "\t" << float_body1->GetAppliedForce().z() << "\n";
			application.DoStep();
			frame++;
		}
		application.EndScene();
	}
	zpos.close();
	auto end = std::chrono::high_resolution_clock::now();
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Duration: " << duration/1000.0 << " seconds" << std::endl;
	return 0;
}
