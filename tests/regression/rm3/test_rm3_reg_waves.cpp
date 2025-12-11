#include <hydroc/gui/guihelper.h>
#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>
#include <hydroc/logging.h>

#include <chrono/core/ChRealtimeStep.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <chrono>   // std::chrono::high_resolution_clock::now
#include <iomanip>  // std::setprecision
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;

int main(int argc, char* argv[]) {
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    // Initialize logging with command line arguments
    // Optional: initialize logging here if test needs debug logs

    // Parse CLI arguments and initialize environment
    bool profilingOn     = true;
    bool saveDataOn      = true;
    bool plotOn          = true;
    bool visualizationOn = true;
    std::string data_dir;
    if (!hydroc::GetCLIArguments(argc, argv, "RM3 regular waves regression test", saveDataOn, profilingOn, plotOn,
                                 visualizationOn, data_dir))
        return 1;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    // Get model file names - use HydroChrono data directory
    std::filesystem::path DATADIR(hydroc::getDataDir());

    std::filesystem::path body1_path = DATADIR / "demos" / "rm3" / "geometry" / "float_cog.obj";
    std::filesystem::path body2_path = DATADIR / "demos" / "rm3" / "geometry" / "plate_cog.obj";
    std::filesystem::path h5_path    = DATADIR / "demos" / "rm3" / "hydroData" / "rm3.h5";

    // Early sanity check so we fail with a clear message if assets are missing
    if (!std::filesystem::exists(body1_path)) {
        std::cerr << "ERROR: RM3 regular-waves test: float_cog mesh not found at " << body1_path << std::endl;
        return 1;
    }
    if (!std::filesystem::exists(body2_path)) {
        std::cerr << "ERROR: RM3 regular-waves test: plate_cog mesh not found at " << body2_path << std::endl;
        return 1;
    }
    if (!std::filesystem::exists(h5_path)) {
        std::cerr << "ERROR: RM3 regular-waves test: rm3.h5 not found at " << h5_path << std::endl;
        return 1;
    }

    auto body1_meshfame = body1_path.lexically_normal().generic_string();
    auto body2_meshfame = body2_path.lexically_normal().generic_string();
    auto h5fname        = h5_path.lexically_normal().generic_string();

    std::cout << "Looking for mesh files in:" << std::endl;
    std::cout << "  body1: " << body1_meshfame << std::endl;
    std::cout << "  body2: " << body2_meshfame << std::endl;
    std::cout << "  h5: " << h5fname << std::endl;

    // system/solver settings
    ChSystemNSC system;

    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
    double timestep = 0.01;
    system.SetTimestepperType(ChTimestepper::Type::HHT);
    system.SetSolverType(ChSolver::Type::GMRES);
    system.GetSolver()->AsIterative()->SetMaxIterations(
        300);  // the higher, the easier to keep the constraints satisfied.
    ChRealtimeStepTimer realtime_timer;
    double simulationDuration = 40.0;

    // some io/viz options
    std::vector<double> time_vector;
    std::vector<double> float_heave_position;
    std::vector<double> float_drift_position;
    std::vector<double> plate_heave_position;

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> float_body1 = chrono_types::make_shared<ChBodyEasyMesh>(  //
        body1_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );
    std::cout << "Created float_body1 from mesh." << std::endl;

    // define the float's initial conditions
    system.Add(float_body1);
    float_body1->SetName("body1");
    float_body1->SetPos(ChVector3d(0, 0, -0.72));
    float_body1->SetMass(725834);
    float_body1->SetInertiaXX(ChVector3d(20907301.0, 21306090.66, 37085481.11));
    // float_body1->EnableCollision(false);

    // Create a visualization material
    auto red = chrono_types::make_shared<ChVisualMaterial>();
    red->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.1f));
    float_body1->GetVisualShape(0)->SetMaterial(0, red);

    std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
    std::shared_ptr<ChBody> plate_body2 = chrono_types::make_shared<ChBodyEasyMesh>(  //
        body2_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );
    std::cout << "Created plate_body2 from mesh." << std::endl;

    // Create a visualization material
    auto blue = chrono_types::make_shared<ChVisualMaterial>();
    blue->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.6f));
    plate_body2->GetVisualShape(0)->SetMaterial(0, blue);

    // define the plate's initial conditions
    system.Add(plate_body2);
    plate_body2->SetName("body2");
    plate_body2->SetPos(ChVector3d(0, 0, (-21.29)));
    plate_body2->SetMass(886691);
    plate_body2->SetInertiaXX(ChVector3d(94419614.57, 94407091.24, 28542224.82));
    // plate_body2->EnableCollision(false);

    // add prismatic joint between the two bodies
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(float_body1, plate_body2, false, ChFramed(ChVector3d(0, 0, -0.72)),
                          ChFramed(ChVector3d(0, 0, -21.29)));
    system.AddLink(prismatic);

    auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
    prismatic_pto->Initialize(float_body1, plate_body2, false, ChVector3d(0, 0, -0.72), ChVector3d(0, 0, -21.29));
    prismatic_pto->SetDampingCoefficient(0.0);
    system.AddLink(prismatic_pto);
    std::cout << "Added prismatic and PTO links." << std::endl;

    // attach hydrodynamic forces to body
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(float_body1);
    bodies.push_back(plate_body2);

    // define wave parameters
    auto my_hydro_inputs                     = std::make_shared<RegularWave>(static_cast<unsigned int>(bodies.size()));
    my_hydro_inputs->regular_wave_amplitude_ = 1.0;
    my_hydro_inputs->regular_wave_omega_     = 2.10;
    std::cout << "Configured RegularWave inputs." << std::endl;

    HydroForces hydro_forces(bodies, h5fname);
    hydro_forces.AddWaves(my_hydro_inputs);

    try {
        hydro_forces.AddWaves(my_hydro_inputs);
        std::cout << "Waves added to TestHydro." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception in TestHydro::AddWaves: " << e.what() << std::endl;
        return 1;
    }
    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "Entering simulation loop..." << std::endl;
    // main simulation loop (headless, no GUI)
    int step = 0;
    while (system.GetChTime() <= simulationDuration) {
        if (step % 100 == 0) {
            std::cout << "Step " << step << ", t = " << system.GetChTime() << std::endl;
        }

        system.DoStepDynamics(timestep);
        ++step;

        // append data to output vector
        time_vector.push_back(system.GetChTime());
        float_heave_position.push_back(float_body1->GetPos().z());
        float_drift_position.push_back(float_body1->GetPos().x());
        plate_heave_position.push_back(plate_body2->GetPos().z());
    }
    std::cout << "Exited simulation loop at t = " << system.GetChTime() << std::endl;

    // for profiling
    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::string out_dir = hydroc::getTestOutDir();
    if (profilingOn || saveDataOn) {
        out_dir = out_dir + "/" + RESULTS_DIR_NAME;
        std::filesystem::create_directory(std::filesystem::path(out_dir));
    }

    if (profilingOn) {
        std::ofstream profilingFile(out_dir + "/" + RESULTS_FILE_NAME + "_duration.txt");
        if (profilingFile.is_open()) {
            profilingFile << duration << " ms\n";
            profilingFile.close();
        } else {
            std::cout << "Error: Could not open profiling file for writing." << std::endl;
        }
    }

    if (saveDataOn) {
        std::ofstream outputFile(out_dir + "/" + RESULTS_FILE_NAME + ".txt");
        if (outputFile.is_open()) {
            outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(16) << "Float Heave (m)"
                       << std::right << std::setw(16) << "Plate Heave (m)" << std::right << std::setw(16)
                       << "Float Drift (x) (m)" << std::endl;
            for (int i = 0; i < time_vector.size(); ++i)
                outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
                           << std::right << std::setw(16) << std::setprecision(4) << std::fixed
                           << float_heave_position[i] << std::right << std::setw(16) << std::setprecision(4)
                           << std::fixed << plate_heave_position[i] << std::right << std::setw(16)
                           << std::setprecision(4) << std::fixed << float_drift_position[i] << std::endl;
            outputFile.close();
        } else {
            std::cout << "Error: Could not open output file for writing." << std::endl;
            return 1;  // Return an error code
        }
    }

    return 0;
}