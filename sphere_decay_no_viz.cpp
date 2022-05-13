#include "hydro_forces.h"
#include "chrono_irrlicht/ChIrrNodeAsset.h"
#include <chrono>

int main(int argc, char* argv[]) {
	auto start = std::chrono::high_resolution_clock::now();
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0, 0, -9.81));

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

	HydroInputs my_hydro_inputs;
	my_hydro_inputs.SetRegularWaveAmplitude(0.022);
	my_hydro_inputs.SetRegularWaveOmega(2.10);
	LoadAllHydroForces blah(body, "../../HydroChrono/sphere.h5", "body1", my_hydro_inputs);

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
		/*if (buttonPressed)*/if(true) {
			zpos << system.GetChTime() << "\t" << body->GetPos().x() << "\t" << body->GetPos().z() << "\t" << body->GetPos_dt().z() << "\t" << body->GetAppliedForce().z() << "\n";
			system.DoStepDynamics(timestep);
			frame++;
		}
	}
	zpos.close();
	auto end = std::chrono::high_resolution_clock::now();
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Duration: " << duration/1000.0 << " seconds" << std::endl;
	return 0;
}
