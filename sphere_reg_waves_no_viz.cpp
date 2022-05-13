#include "hydro_forces.h"
#include "chrono_irrlicht/ChIrrNodeAsset.h"
#include <filesystem>
#include <chrono>

int main(int argc, char* argv[]) {
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";
	GetLog() << "HydroChrono v0.0.1\n\n";

	int regWaveNum = 10;

	double task10WaveAmps[] = { 0.044, 0.078, 0.095, 0.123, 0.177, 0.24, 0.314, 0.397, 0.491, 0.594 };
	double task10WaveOmegas[] = { 2.094395102, 1.570796327, 1.427996661, 1.256637061, 1.047197551, 0.897597901, 0.785398163, 0.698131701, 0.628318531, 0.571198664 };
	double task10dampings[] = { 398736.034, 118149.758, 90080.857, 161048.558, 322292.419, 479668.979, 633979.761, 784083.286, 932117.647, 1077123.445 };
	//int waveNum = 0;

	//std::cout << '\n' << "Enter regWaveNum (1-10)...\n";
	//std::cin >> regWaveNum;

	// define some basic model parameters
	auto start = std::chrono::high_resolution_clock::now();
	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0, 0, -9.81));

	// Setup Ground
	auto ground = chrono_types::make_shared<ChBody>();
	system.AddBody(ground);
	ground->SetPos(ChVector<>(0, 0, -5));
	ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(false);

	// create easy sphere
	std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasySphere>(5, 1);
	auto sph = chrono_types::make_shared<ChSphereShape>();
	body->AddAsset(sph);
	system.Add(body);
	// set up body initial conditions
	body->SetPos(ChVector<>(0, 0, -2));
	body->SetMass(261.8e3);
	// attach color asset to body
	auto col_2 = chrono_types::make_shared<ChColorAsset>();
	col_2->SetColor(ChColor(0, 0, 0.6f));
	body->AddAsset(col_2);

	// set up output file for body position each step
	std::string out_dir = "results/regular_waves/";
	std::cout << "creating results directory...\n";
	std::filesystem::create_directories(out_dir);

	// S = 0.0005
	//std::cout << "running wave no. " + std::to_string(regWaveNum) + "\n";
	// Info about which solver to use - may want to change this later
	auto gmres_solver = chrono_types::make_shared<ChSolverGMRES>();  // change to mkl or minres?
	gmres_solver->SetMaxIterations(300);
	system.SetSolver(gmres_solver);
	double timestep = 0.015; // also sets the timesteps in chrono system

	// Create the spring between body_1 and ground. The spring end points are
	// specified in the body relative frames.
	double rest_length = 3.0;
	double spring_coef = 0.0;
	double damping_coef = task10dampings[regWaveNum - 1];
	auto spring_1 = chrono_types::make_shared<ChLinkTSDA>();
	spring_1->Initialize(body, ground, true, ChVector<>(0, 0, -2), ChVector<>(0, 0, -5));
	spring_1->SetRestLength(rest_length);
	spring_1->SetSpringCoefficient(spring_coef);
	spring_1->SetDampingCoefficient(damping_coef);
	system.AddLink(spring_1);
	//// attach color asset to spring
	//auto col_1 = chrono_types::make_shared<ChColorAsset>();
	//col_1->SetColor(ChColor(0, 0, 0));
	//spring_1->AddAsset(col_1);
	//// Attach a visualization asset.
	//spring_1->AddAsset(col_1);
	//spring_1->AddAsset(chrono_types::make_shared<ChPointPointSpring>(2, 80, 15));

	HydroInputs myHydroInputs;
	myHydroInputs.regularWaveAmplitude = task10WaveAmps[regWaveNum-1]; //0.095;
	myHydroInputs.regularWaveOmega = task10WaveOmegas[regWaveNum-1];//1.427996661;
	LoadAllHydroForces blah(body, "../../HydroChrono/sphere.h5", "body1", myHydroInputs);

	std::string out_file = "regwave_" + std::to_string(regWaveNum) + ".txt";
	std::ofstream out_stream(out_dir + out_file, std::ofstream::out);
	//if (!out_stream.is_open()) {
	//	std::cout << "Error opening file \"" + out_file + "\". Please make sure this file path exists then try again\n";
	//	return -1;
	//}
	out_stream.precision(10);
	out_stream.width(12);
	out_stream << "Wave #: \t" << regWaveNum << "\n";
	out_stream << "Wave amplitude (m): \t" << myHydroInputs.regularWaveAmplitude << "\n";
	out_stream << "Wave omega (rad/s): \t" << myHydroInputs.regularWaveOmega << "\n";
	out_stream << "#Time\tBody Pos\tBody vel (heave)\tforce (heave)\n";

	// Simulation loop
	int frame = 0;
	while (system.GetChTime() <= 400) {
		if (true) {
			out_stream << system.GetChTime() << "\t" << body->GetPos().x() << "\t" << body->GetPos().z() << "\t" << body->GetPos_dt().z() << "\t" << body->GetAppliedForce().z() << "\n";
			system.DoStepDynamics(timestep);
			frame++;
		}
	}
	out_stream.close();
	auto end = std::chrono::high_resolution_clock::now();
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Duration: " << duration / 1000.0 << " seconds" << std::endl;

	//do
	//{
	//	std::cout << '\n' << "Press any key to continue...";
	//} while (std::cin.get() != '\n');
	return 0;
}