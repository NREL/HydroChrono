#include <hydroc/hydro_forces.h>
#include <hydroc/helper.h>
#include <hydroc/gui/guihelper.h>

#include <chrono/core/ChRealtimeStep.h>

#include <filesystem>
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

	auto body1_meshfame = (DATADIR / "sphere" / "geometry" /"oes_task10_sphere.obj")
		.lexically_normal().generic_string();
	auto h5fname = (DATADIR / "sphere" / "hydroData" /"sphere.h5")
		.lexically_normal().generic_string();



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
	std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);

	hydroc::gui::UI& ui = *pui.get();

	// Setup Ground
	auto ground = chrono_types::make_shared<ChBody>();
	system.AddBody(ground);
	ground->SetPos(ChVector<>(0, 0, -5));
	ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(false);

	// some io/viz options
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
	system.Add(sphereBody);
	sphereBody->SetNameString("body1"); // must set body name correctly! (must match .h5 file)
	sphereBody->SetPos(ChVector<>(0, 0, -2));
	sphereBody->SetMass(261.8e3);

	// Create a visualization material
	auto red = chrono_types::make_shared<ChVisualMaterial>();
	red->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.1f));
	sphereBody->GetVisualShape(0)->SetMaterial(0, red);


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
	my_hydro_inputs.mode = WaveMode::regular; // uses regular wave mode
	my_hydro_inputs.regular_wave_amplitude = task10_wave_amps[reg_wave_num-1]; //0.095;
	my_hydro_inputs.regular_wave_omega = task10_wave_omegas[reg_wave_num-1];//1.427996661;

	std::vector<std::shared_ptr<ChBody>> bodies;
	bodies.push_back(sphereBody);
	TestHydro blah(bodies, h5fname, my_hydro_inputs);

	// for profiling
	auto start = std::chrono::high_resolution_clock::now();
	
	// main simulation loop
	ui.Init(&system, "Sphere Reg Waves - Decay Test"); 
	ui.SetCamera(8, -25, 15, 0, 0, 0);

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
		std::string out_dir = "results/sphere_regular_waves/";
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
		std::string out_dir = "results/sphere_regular_waves/";
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
