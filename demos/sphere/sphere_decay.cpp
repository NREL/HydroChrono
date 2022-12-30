#include <hydroc/hydro_forces.h>
#include <hydroc/helper.h>

#include <hydroc/gui/guihelper.h>


#include "chrono/core/ChRealtimeStep.h"


#include <iomanip> // std::setprecision
#include <chrono> // std::chrono::high_resolution_clock::now
#include <vector> // std::vector<double>
#include <filesystem> // c++17 only

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;



// the main program to be executed:
int main(int argc, char* argv[]) {
	GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";

    if (hydroc::setInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    std::filesystem::path DATADIR(hydroc::getDataDir());

	auto body1_meshfame = (DATADIR / "sphere" / "geometry" /"oes_task10_sphere.obj").lexically_normal().generic_string();
	auto h5fname = (DATADIR / "sphere" / "hydroData" /"sphere.h5").lexically_normal().generic_string();



	// system/solver settings
	ChSystemNSC system;

	system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));

	double timestep = 0.015;
	system.SetSolverType(ChSolver::Type::GMRES);
	system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
	system.SetStep(timestep);
	ChRealtimeStepTimer realtime_timer;
	double simulationDuration = 40.0;

	// Create user interface
	bool visualizationOn = true;

	std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);

	hydroc::gui::UI& ui = *pui.get();

	bool profilingOn = true;
	bool saveDataOn = true;

	// Output timeseries
	std::vector<double> time_vector;
	std::vector<double> heave_position;

	// set up body from a mesh
	std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
	std::shared_ptr<ChBody> sphereBody = chrono_types::make_shared<ChBodyEasyMesh>(       //
		body1_meshfame,    // file name
		1000,                                                                             // density
		false,                                                                            // do not evaluate mass automatically
		true,                                                                             // create visualization asset
		false                                                                             // do not collide
		);

	// define the body's initial conditions
	sphereBody->SetNameString("body1"); // must set body name correctly! (must match .h5 file)
	sphereBody->SetPos(ChVector<>(0, 0, -1));
	sphereBody->SetMass(261.8e3);

	// Create a visualization material
	auto cadet_blue = chrono_types::make_shared<ChVisualMaterial>();
	cadet_blue->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.1f));
	sphereBody->GetVisualShape(0)->SetMaterial(0, cadet_blue);

	system.Add(sphereBody);


	// define wave parameters (not used in this demo)
	// Todo define a way to use TestHydro without hydro_inputs/waves
	HydroInputs my_hydro_inputs;
	my_hydro_inputs.mode = WaveMode::noWaveCIC;
	//my_hydro_inputs.regular_wave_amplitude = 0.022;
	//my_hydro_inputs.regular_wave_omega = 2.10;

	// attach hydrodynamic forces to body
	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(sphereBody);


	TestHydro blah(bodies, h5fname, my_hydro_inputs);

	// for profilingvisualizationOn = false;
	auto start = std::chrono::high_resolution_clock::now(); 

	// main simulation loop
	ui.Init(&system, "Sphere - Decay Test"); 

	while (system.GetChTime() <= simulationDuration) {

		if(ui.IsRunning(timestep) == false) break;
		
		if (ui.simulationStarted) {

			// append data to output vector
			time_vector.push_back(system.GetChTime());
			heave_position.push_back(sphereBody->GetPos().z());

		}
	}


	// for profiling
	auto end = std::chrono::high_resolution_clock::now();
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	if (profilingOn) {
		std::ofstream profilingFile;
		profilingFile.open("./results/decay/duration_ms.txt");
		if (!profilingFile.is_open()) {
			if (!std::filesystem::exists("./results/decay")) {
				std::cout << "Path " << std::filesystem::absolute("./results/decay") << " does not exist, creating it now..." << std::endl;
				std::filesystem::create_directory("./results");
				std::filesystem::create_directory("./results/decay");
				profilingFile.open("./results/decay/duration_ms.txt");
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
		outputFile.open("./results/decay/sphere_decay.txt");
		if (!outputFile.is_open()) {
			if (!std::filesystem::exists("./results/decay")) {
				std::cout << "Path " << std::filesystem::absolute("./results/decay") << " does not exist, creating it now..." << std::endl;
				std::filesystem::create_directory("./results");
				std::filesystem::create_directory("./results/decay");
				outputFile.open("./results/decay/sphere_decay.txt");
				if (!outputFile.is_open()) {
					std::cout << "Still cannot open file, ending program" << std::endl;
					return 0;
				}
			}
		}
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

	std::cout << "Simulation finished." << std::endl;
	return 0;
}