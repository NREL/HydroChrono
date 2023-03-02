#include <hydroc/hydro_forces.h>

#include <hydroc/helper.h>

#ifdef HYDROCHRONO_HAVE_IRRLICHT
    #include <chrono_irrlicht/ChVisualSystemIrrlicht.h>
    #include <chrono_irrlicht/ChIrrMeshTools.h>
    // Use the main namespaces of Irrlicht
    using namespace irr;
    using namespace irr::core;
    using namespace irr::scene;
    using namespace irr::video;
    using namespace irr::io;
    using namespace irr::gui;
    using namespace chrono::irrlicht;
#endif

#include <chrono/core/ChRealtimeStep.h>
#include <chrono/physics/ChLinkMate.h>
#include <iomanip> // std::setprecision
#include <chrono> // std::chrono::high_resolution_clock::now
#include <vector> // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;


#ifdef HYDROCHRONO_HAVE_IRRLICHT
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
#endif


int main(int argc, char* argv[]) {
	//auto start = std::chrono::high_resolution_clock::now();
	GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";


	if (hydroc::setInitialEnvironment(argc, argv) != 0) {
		return 1;
	}

	std::filesystem::path DATADIR(hydroc::getDataDir());

	auto body1_meshfame = (DATADIR / "DeepCWind" / "geometry" / "deepcwind.obj")
		.lexically_normal()
		.generic_string();
	auto h5fname = (DATADIR / "DeepCWind" / "hydroData" / "deepcwind.h5")
		.lexically_normal()
		.generic_string();

	// system/solver settings
	ChSystemSMC system;
	system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
	double timestep = 0.08;
	system.SetTimestepperType(ChTimestepper::Type::HHT);
	system.SetSolverType(ChSolver::Type::GMRES);
	system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
	system.SetStep(timestep);
	ChRealtimeStepTimer realtime_timer;
	double simulationDuration = 1000.0;

	// some io/viz options
	bool visualizationOn = true;
	bool profilingOn     = true;
	bool saveDataOn      = true;
	std::vector<double> time_vector;
	std::vector<double> base_pitch;
	std::vector<double> base_surge;


	// set up base from a mesh
	std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
	std::shared_ptr<ChBody> base = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		body1_meshfame,
		0,                                                                                        // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false                                                                                     // collisions
		);

	// define the base's initial conditions 
	system.Add(base);
	base->SetNameString("body1");
	auto cg     = ChVector<>(0.0, 0.0, -7.53);
	// offset used for heave/surge decay test
	auto offset = ChVector<>(0.0, 0.0, 0.0);
	base->SetPos(cg + offset);
	// Use for pitch decay test
	double ang_rad = -3.95 * CH_C_PI / 180.0;
	base->SetRot(Q_from_AngAxis(ang_rad, VECT_Y));
	base->SetMass(1.419625e7);
	base->SetInertiaXX(ChVector<>(1.2898e10, 1.2851e10, 1.4189e10)); 

	// add fixed ground for linear damping (surge or pitch)
	auto ground = chrono_types::make_shared<ChBody>();
	system.AddBody(ground);
	ground->SetPos(cg);
	ground->SetRot(Q_from_AngAxis(CH_C_PI / 2.0, VECT_X));
	ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(false);

	// define damping in pitch
	auto rot_damp = chrono_types::make_shared<ChLinkRSDA>();
	// need to set damping to 31 MN-m/(rad/s)
	rot_damp->SetDampingCoefficient(31e6);
	// puts Z axis for link into screen, keeping x axis the same (to the right)
	ChQuaternion<> rev_rot = Q_from_AngAxis(CH_C_PI / 2.0, VECT_X); // do not change
	rot_damp->Initialize(base, ground, false, ChCoordsys(cg, rev_rot), ChCoordsys(cg, rev_rot));
	system.AddLink(rot_damp);

	// define wave parameters (not used in this demo TODO have hydroforces constructor without hydro inputs)
	HydroInputs my_hydro_inputs;
	my_hydro_inputs.mode = WaveMode::noWaveCIC;
	
	// set up hydro forces
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(base);
	TestHydro hydroforces(bodies, h5fname, my_hydro_inputs);

	// for profiling
	auto start = std::chrono::high_resolution_clock::now();

#ifdef HYDROCHRONO_HAVE_IRRLICHT
	if (visualizationOn) {
		// create the irrlicht application for visualizing
		auto irrlichtVis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
		irrlichtVis->AttachSystem(&system);
		irrlichtVis->SetWindowSize(1280, 720);
		irrlichtVis->SetWindowTitle("DeepCwind Verification");
		irrlichtVis->SetCameraVertical(CameraVerticalDir::Z);
		irrlichtVis->Initialize();
		irrlichtVis->AddLogo();
		irrlichtVis->AddSkyBox();
		irrlichtVis->AddCamera(ChVector<>(0, -70, -10), ChVector<>(0, 0, -10));
		irrlichtVis->AddTypicalLights();
		//irrlichtVis->EnableBodyFrameDrawing(true);
		//irrlichtVis->EnableLinkFrameDrawing(true);

		// add play/pause button
		bool buttonPressed = false;
		MyActionReceiver receiver(irrlichtVis.get(), buttonPressed);
		irrlichtVis->AddUserEventReceiver(&receiver);
		//ChSparseMatrix M;
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
				base_surge.push_back(base->GetPos().x());
				base_pitch.push_back(base->GetRot().Q_to_Euler123().y());
				// force playback to be real-time
				realtime_timer.Spin(timestep);

			}
		}
	}
	else {
#endif // #ifdef HYDROCHRONO_HAVE_IRRLICHT
		while (system.GetChTime() <= simulationDuration) {
			// append data to std vector
			time_vector.push_back(system.GetChTime());
			base_surge.push_back(base->GetPos().x());
			base_pitch.push_back(base->GetRot().Q_to_Euler123().y());
			// step the simulation forwards
			system.DoStepDynamics(timestep);
		}
#ifdef HYDROCHRONO_HAVE_IRRLICHT
	}
#endif

	if (saveDataOn) {
		std::ofstream outputFile;
		outputFile.open("./results/DeepCWind/decay/DeepCWind_decay.txt");
		if (!outputFile.is_open()) {
			if (!std::filesystem::exists("./results/DeepCWind/decay")) {
				std::cout << "Path " << std::filesystem::absolute("./results/DeepCWind/decay") << " does not exist, creating it now..." << std::endl;
				std::filesystem::create_directory("./results");
				std::filesystem::create_directory("./results/DeepCWind");
				std::filesystem::create_directory("./results/DeepCWind/decay");
				outputFile.open("./results/DeepCWind/decay/DeepCWind_decay.txt");
				if (!outputFile.is_open()) {
					std::cout << "Still cannot open file, ending program" << std::endl;
					return 0;
				}
			}
		}
		outputFile << std::left << std::setw(10) << "Time (s)"
			<< std::right << std::setw(16) << "Base Surge (m)"
			<< std::right << std::setw(16) << "Base Pitch (radians)"
			<< std::endl;
		for (int i = 0; i < time_vector.size(); ++i)
			outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
			<< std::right << std::setw(16) << std::setprecision(8) << std::fixed << base_surge[i]
			<< std::right << std::setw(16) << std::setprecision(8) << std::fixed << base_pitch[i]
			<< std::endl;
		outputFile.close();
	}
	return 0;
}