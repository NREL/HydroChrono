#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <chrono>   // std::chrono::high_resolution_clock::now
#include <iomanip>  // std::setprecision
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;

int main(int argc, char* argv[]) {
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    // Initialize environment
    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    // Get model file names
    std::filesystem::path DATADIR(hydroc::getDataDir());

    std::filesystem::path body1_path = DATADIR / "demos" / "rm3" / "geometry" / "float_cog.obj";
    std::filesystem::path body2_path = DATADIR / "demos" / "rm3" / "geometry" / "plate_cog.obj";
    std::filesystem::path h5_path    = DATADIR / "demos" / "rm3" / "hydroData" / "rm3.h5";

    // Early sanity check
    if (!std::filesystem::exists(body1_path)) {
        std::cerr << "ERROR: float_cog mesh not found at " << body1_path << std::endl;
        return 1;
    }
    if (!std::filesystem::exists(body2_path)) {
        std::cerr << "ERROR: plate_cog mesh not found at " << body2_path << std::endl;
        return 1;
    }
    if (!std::filesystem::exists(h5_path)) {
        std::cerr << "ERROR: rm3.h5 not found at " << h5_path << std::endl;
        return 1;
    }

    auto body1_meshfame = body1_path.lexically_normal().generic_string();
    auto body2_meshfame = body2_path.lexically_normal().generic_string();
    auto h5fname        = h5_path.lexically_normal().generic_string();

    // system/solver settings
    ChSystemNSC system;

    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));

    double timestep = 0.01;
    system.SetTimestepperType(ChTimestepper::Type::HHT);
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
    double simulationDuration = 40.0;

    // Output timeseries
    std::vector<double> time_vector;
    std::vector<double> float_heave_position;
    std::vector<double> plate_heave_position;

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> float_body1 = chrono_types::make_shared<ChBodyEasyMesh>(
        body1_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the float's initial conditions
    system.Add(float_body1);
    float_body1->SetName("body1");
    float_body1->SetPos(ChVector3d(0, 0, (-0.72 + 0.1)));
    float_body1->SetMass(725834);
    float_body1->SetInertiaXX(ChVector3d(20907301.0, 21306090.66, 37085481.11));

    // Plate
    std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
    std::shared_ptr<ChBody> plate_body2 = chrono_types::make_shared<ChBodyEasyMesh>(
        body2_meshfame,
        0,      // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the plate's initial conditions
    system.Add(plate_body2);
    plate_body2->SetName("body2");
    plate_body2->SetPos(ChVector3d(0, 0, (-21.29)));
    plate_body2->SetMass(886691);
    plate_body2->SetInertiaXX(ChVector3d(94419614.57, 94407091.24, 28542224.82));

    // add prismatic joint between the two bodies
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(float_body1, plate_body2, false, ChFramed(ChVector3d(0, 0, -0.72)),
                          ChFramed(ChVector3d(0, 0, -21.29)));
    system.AddLink(prismatic);

    auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
    prismatic_pto->Initialize(float_body1, plate_body2, false, ChVector3d(0, 0, -0.72), ChVector3d(0, 0, -21.29));
    prismatic_pto->SetDampingCoefficient(0.0);
    system.AddLink(prismatic_pto);

    auto default_dont_add_waves = std::make_shared<NoWave>(2);

    // attach hydrodynamic forces to body
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(float_body1);
    bodies.push_back(plate_body2);

    HydroForces hydro_forces(bodies, h5fname);
    hydro_forces.AddWaves(default_dont_add_waves);

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

    // main simulation loop
    while (system.GetChTime() <= simulationDuration) {
        system.DoStepDynamics(timestep);

        // append data to output vector
        time_vector.push_back(system.GetChTime());
        float_heave_position.push_back(float_body1->GetPos().z());
        plate_heave_position.push_back(plate_body2->GetPos().z());
    }

    // for profiling
    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Simulation completed in " << duration << " ms" << std::endl;

    // Save results
    std::string out_dir = hydroc::getTestOutDir() + "/" + RESULTS_DIR_NAME;
    std::filesystem::create_directories(std::filesystem::path(out_dir));

    std::ofstream outputFile(out_dir + "/" + RESULTS_FILE_NAME + ".txt");
    if (outputFile.is_open()) {
        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(16) << "Float Heave (m)"
                   << std::right << std::setw(16) << "Plate Heave (m)" << std::endl;
        for (size_t i = 0; i < time_vector.size(); ++i)
            outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
                       << std::right << std::setw(16) << std::setprecision(8) << std::fixed
                       << float_heave_position[i] << std::right << std::setw(16) << std::setprecision(8)
                       << std::fixed << plate_heave_position[i] << std::endl;
        outputFile.close();
    } else {
        std::cout << "Error: Could not open output file for writing." << std::endl;
        return 1;
    }

    return 0;
}
