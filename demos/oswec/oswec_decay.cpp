#include <hydroc/hydro_forces.h>
#include <hydroc/helper.h>
#include <hydroc/gui/guihelper.h>

#include <chrono/core/ChRealtimeStep.h>

#include <iomanip> // std::setprecision
#include <chrono> // std::chrono::high_resolution_clock::now
#include <vector> // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;

// usage: ./<demos>.exe [DATADIR] [--nogui]
//
// If no argument is given user can set HYDROCHRONO_DATA_DIR 
// environment variable to give the data_directory.
// 
int main(int argc, char* argv[]) {

	GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";

	if (hydroc::setInitialEnvironment(argc, argv) != 0) {
		return 1;
	}

	// Check if --nogui option is set as 2nd argument
	bool visualizationOn = true;
	if(argc > 2 &&  std::string("--nogui").compare(argv[2]) == 0)  {
		visualizationOn = false;
	}

	// Get model file names
	std::filesystem::path DATADIR(hydroc::getDataDir());

	auto body1_meshfame = (DATADIR / "oswec" / "geometry" / "flap.obj")
		.lexically_normal()
		.generic_string();
	auto body2_meshfame = (DATADIR / "oswec" / "geometry" / "base.obj")
		.lexically_normal()
		.generic_string();
	auto h5fname = (DATADIR / "oswec" / "hydroData" / "oswec.h5")
		.lexically_normal()
		.generic_string();

	// system/solver settings
	ChSystemNSC system;

	system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
	double timestep = 0.03;
	//system.SetTimestepperType(ChTimestepper::Type::HHT);
	system.SetSolverType(ChSolver::Type::GMRES);
	//system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
	system.SetStep(timestep);
	ChRealtimeStepTimer realtime_timer;
	double simulationDuration = 400.0;

	// Create user interface
	std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);
	hydroc::gui::UI& ui = *pui.get();

	// some io/viz options
	bool profilingOn = false;
	bool saveDataOn = true;
	std::vector<double> time_vector;
	std::vector<double> flap_rot;

	// set up body from a mesh
	std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
	std::shared_ptr<ChBody> flap_body = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		body1_meshfame,
		1000,                                                                                     // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false                                                                                     // collisions
		);

	// Create a visualization material
	auto red = chrono_types::make_shared<ChVisualMaterial>();
	red->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.1f));
	flap_body->GetVisualShape(0)->SetMaterial(0, red);

	// define the float's initial conditions
	system.Add(flap_body);
	flap_body->SetNameString("body1");
	auto ang_rad = CH_C_PI / 18.0;
	flap_body->SetPos(ChVector<>(6.1*std::cos(CH_C_PI / 2.0 - ang_rad),
		                       0.0,
		                      -10.0+6.1*std::sin(CH_C_PI / 2.0 -ang_rad)));
	flap_body->SetRot(Q_from_AngAxis(ang_rad, VECT_Y));
	flap_body->SetMass(127000.0);
	flap_body->SetInertiaXX(ChVector<>(1.85e6, 1.85e6, 1.85e6));
	// notes: mass and inertia added to added mass and system mass correctly.


	// set up body from a mesh
	std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
	std::shared_ptr<ChBody> base_body = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		body2_meshfame,
		1000,                                                                                     // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false                                                                                     // collisions
		);

	// Create a visualization material
	auto blue = chrono_types::make_shared<ChVisualMaterial>();
	blue->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.6f));
	base_body->GetVisualShape(0)->SetMaterial(0, blue);

	// define the plate's initial conditions
	system.Add(base_body);
	base_body->SetNameString("body2");
	base_body->SetPos(ChVector<>(0, 0, -10.9));
	base_body->SetMass(999);
	base_body->SetInertiaXX(ChVector<>(1, 1, 1));
	base_body->SetBodyFixed(true);

	// define base-fore flap joint
	ChQuaternion<> revoluteRot = Q_from_AngX(CH_C_PI / 2.0);
	auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
	revolute->Initialize(base_body, flap_body, ChCoordsys<>(ChVector<>(0.0, 0.0, -10.0), revoluteRot));
	system.AddLink(revolute);

	//// define wave parameters (not used in this demo)
	HydroInputs my_hydro_inputs;
	my_hydro_inputs.mode = WaveMode::noWaveCIC;// or 'regular' or 'regularCIC' or 'irregular';
	////my_hydro_inputs.regular_wave_amplitude = 0.022;
	////my_hydro_inputs.regular_wave_omega = 2.10;

	//// attach hydrodynamic forces to body
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(flap_body);
	bodies.push_back(base_body);
	TestHydro blah(bodies, h5fname, my_hydro_inputs);

	// for profiling
	auto start = std::chrono::high_resolution_clock::now();

	// main simulation loop
	ui.Init(&system, "OSWEC - Decay Test");
	ui.SetCamera(0, -50, -10, 0, 0, -10);

	while (system.GetChTime() <= simulationDuration) {

		if(ui.IsRunning(timestep) == false) break;
		
		if (ui.simulationStarted) {

			// append data to output vector
			time_vector.push_back(system.GetChTime());
			flap_rot.push_back(flap_body->GetRot().Q_to_Euler123().y());		
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
			<< std::right << std::setw(16) << "Flap Rotation y (radians)"
			<< std::right << std::setw(16) << "Flap Rotation y (degrees)"
			<< std::endl;
		for (int i = 0; i < time_vector.size(); ++i)
			outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
			<< std::right << std::setw(16) << std::setprecision(4) << std::fixed << flap_rot[i]
			<< std::right << std::setw(16) << std::setprecision(4) << std::fixed << flap_rot[i]*360.0/6.28
			<< std::endl;
		outputFile.close();
	}
	return 0;
}