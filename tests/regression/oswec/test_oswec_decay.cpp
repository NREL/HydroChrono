#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <chrono>   // std::chrono::high_resolution_clock::now
#include <iomanip>  // std::setprecision
#include <vector>   // std::vector<double>

// Use the namespaces of Chrono
using namespace chrono;

// Function to compute cross product
std::array<double, 3> cross(std::array<double, 3> v1, std::array<double, 3> v2) {
    return {v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]};
}

// Function to compute dot product
double dot(std::array<double, 3> v1, std::array<double, 3> v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Function to normalize a vector
std::array<double, 3> normalize(std::array<double, 3> v) {
    double norm = sqrt(dot(v, v));
    return {v[0] / norm, v[1] / norm, v[2] / norm};
}

// Function to rotate a vector
std::array<double, 3> rotate_vector_3d(std::array<double, 3> vector,
                                       std::array<double, 3> axis,
                                       double angle_in_degrees) {
    double angle_in_radians = angle_in_degrees * CH_DEG_TO_RAD;
    axis = normalize(axis);
    std::array<double, 3> rotated_vector;
    for (int i = 0; i < 3; i++) {
        rotated_vector[i] = vector[i] * cos(angle_in_radians) + cross(axis, vector)[i] * sin(angle_in_radians) +
                            axis[i] * dot(axis, vector) * (1 - cos(angle_in_radians));
    }
    return rotated_vector;
}

// Function to add two vectors
std::array<double, 3> add_vectors(std::array<double, 3> v1, std::array<double, 3> v2) {
    return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

int main(int argc, char* argv[]) {
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    // Initialize environment
    std::string data_dir;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    // Get model file names
    std::filesystem::path DATADIR(hydroc::getDataDir());

    auto body1_meshfame = (DATADIR / "demos" / "oswec" / "geometry" / "flap.obj").lexically_normal().generic_string();
    auto body2_meshfame = (DATADIR / "demos" / "oswec" / "geometry" / "base.obj").lexically_normal().generic_string();
    auto h5fname        = (DATADIR / "demos" / "oswec" / "hydroData" / "oswec.h5").lexically_normal().generic_string();

    // system/solver settings
    ChSystemNSC system;

    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
    double timestep = 0.03;
    system.SetSolverType(ChSolver::Type::GMRES);
    double simulationDuration = 400.0;

    // Output timeseries
    std::vector<double> time_vector;
    std::vector<double> flap_rot;

    std::array<double, 3> origin_to_hinge = {0, 0, -8.9};
    std::array<double, 3> hinge_to_cg     = {0, 0, 5};
    std::array<double, 3> axis            = {0, 1, 0};
    double angle_in_degrees               = 10;

    std::array<double, 3> rotated_hinge_to_cg = rotate_vector_3d(hinge_to_cg, axis, angle_in_degrees);
    std::array<double, 3> new_cg = add_vectors(origin_to_hinge, rotated_hinge_to_cg);

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> flap_body = chrono_types::make_shared<ChBodyEasyMesh>(
        body1_meshfame,
        1000,   // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the float's initial conditions
    system.Add(flap_body);
    flap_body->SetName("body1");
    auto ang_rad = CH_PI / 18.0;
    flap_body->SetPos(ChVector3d(new_cg[0], new_cg[1], new_cg[2]));
    flap_body->SetRot(QuatFromAngleY(ang_rad));
    flap_body->SetMass(127000.0);
    flap_body->SetInertiaXX(ChVector3d(1.85e6, 1.85e6, 1.85e6));

    // set up body from a mesh
    std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
    std::shared_ptr<ChBody> base_body = chrono_types::make_shared<ChBodyEasyMesh>(
        body2_meshfame,
        1000,   // density
        false,  // do not evaluate mass automatically
        true,   // create visualization asset
        false   // collisions
    );

    // define the plate's initial conditions
    system.Add(base_body);
    base_body->SetName("body2");
    base_body->SetPos(ChVector3d(0, 0, -10.15));
    base_body->SetMass(999);
    base_body->SetInertiaXX(ChVector3d(1, 1, 1));

    // create ground
    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetPos(ChVector3d(0, 0, -10.15));
    ground->SetTag(-1);
    ground->SetFixed(true);
    ground->EnableCollision(false);

    // fix base to ground with special constraint
    auto anchor = chrono_types::make_shared<ChLinkMateGeneric>();
    anchor->Initialize(base_body, ground, false, base_body->GetVisualModelFrame(), base_body->GetVisualModelFrame());
    system.Add(anchor);
    anchor->SetConstrainedCoords(true, true, true, true, true, true);

    // define base-fore flap joint
    ChQuaternion<> revoluteRot = QuatFromAngleX(CH_PI / 2.0);
    auto revolute              = chrono_types::make_shared<ChLinkLockRevolute>();
    revolute->Initialize(base_body, flap_body, ChFramed(ChVector3d(0.0, 0.0, -8.9), revoluteRot));
    system.AddLink(revolute);

    auto default_dont_add_waves = std::make_shared<NoWave>(2);

    // set up hydro forces
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(flap_body);
    bodies.push_back(base_body);
    HydroForces blah(bodies, h5fname, default_dont_add_waves);

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

    // main simulation loop
    while (system.GetChTime() <= simulationDuration) {
        system.DoStepDynamics(timestep);

        // append data to output vector
        time_vector.push_back(system.GetChTime());
        flap_rot.push_back(flap_body->GetRot().GetCardanAnglesXYZ().y());
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
        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(16)
                   << "Flap Rotation y (radians)" << std::right << std::setw(16) << "Flap Rotation y (degrees)"
                   << std::endl;
        for (size_t i = 0; i < time_vector.size(); ++i)
            outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
                       << std::right << std::setw(16) << std::setprecision(4) << std::fixed << flap_rot[i]
                       << std::right << std::setw(16) << std::setprecision(4) << std::fixed
                       << flap_rot[i] * 360.0 / 6.28 << std::endl;
        outputFile.close();
    } else {
        std::cout << "Error: Could not open output file for writing." << std::endl;
        return 1;
    }

    return 0;
}
