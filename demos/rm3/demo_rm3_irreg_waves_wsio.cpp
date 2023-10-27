#include <hydroc/gui/guihelper.h>
#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/core/ChRealtimeStep.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;

struct SimulationConfig {
    std::string explorer;
    double startTime;
    double rampTime;
    double endTime;
    double dt;
};

struct WaveConfig {
    std::string type;
    double height;
    double period;
    std::string spectrumType;
    std::vector<double> direction;
    std::string elevationFile;
};

struct BodyConfig {
    std::string hydroDataFile;
    std::string geometryFile;
    double mass;
    std::vector<double> inertia;
    std::vector<double> position;  // {x, y, z}
};

struct PTOConfig {
    std::string type;
    double stiffness;
    double damping;
    std::vector<double> location;
    std::vector<int> bodies;
    std::vector<std::vector<double>> attachments;
};

// Function to parse the Simulation Configuration
SimulationConfig parseSimulationConfig(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        throw std::runtime_error("Failed to open file.");
    }

    SimulationConfig config;
    std::string line;
    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        std::string token;

        if (line.find("simu.explorer = '") != std::string::npos) {
            config.explorer = line.substr(line.find("'") + 1, line.rfind("'") - line.find("'") - 1);
        } else if (line.find("simu.startTime =") != std::string::npos) {
            iss >> token >> token >> config.startTime;
        } else if (line.find("simu.rampTime =") != std::string::npos) {
            iss >> token >> token >> config.rampTime;
        } else if (line.find("simu.endTime =") != std::string::npos) {
            iss >> token >> token >> config.endTime;
        } else if (line.find("simu.dt =") != std::string::npos) {
            iss >> token >> token >> config.dt;
        }
    }

    std::cout << "Explorer: " << config.explorer << std::endl;
    std::cout << "Start Time: " << config.startTime << std::endl;
    std::cout << "Ramp Time: " << config.rampTime << std::endl;
    std::cout << "End Time: " << config.endTime << std::endl;
    std::cout << "Time Step (dt): " << config.dt << std::endl;

    return config;
}

WaveConfig parseWaveConfig(const std::string& filePath) {
    std::ifstream file(filePath);
    std::string line;
    WaveConfig waveConfig;

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return waveConfig;
    }

    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string word;

        iss >> word;
        if (word == "waves") {
            iss >> word;  // Getting '='
            std::string temp;
            iss >> temp;  // Getting the full string, e.g., waveClass('regular');
            // Extracting the word within single quotes
            size_t start    = temp.find("'") + 1;
            size_t end      = temp.find_last_of("'");
            waveConfig.type = temp.substr(start, end - start);
        } else if (word == "waves.height") {
            iss >> word;  // Getting '='
            iss >> waveConfig.height;
        } else if (word == "waves.period") {
            iss >> word;  // Getting '='
            iss >> waveConfig.period;
        } else if (word == "waves.spectrumType") {
            iss >> word;  // Getting '='
            iss >> waveConfig.spectrumType;
        } else if (word == "waves.direction") {
            iss >> word;  // Getting '='
            double direction;
            while (iss >> direction) {
                waveConfig.direction.push_back(direction);
            }
        } else if (word == "waves.elevationFile") {
            iss >> word;  // Getting '='
            iss >> waveConfig.elevationFile;
        }
    }

    file.close();

    std::cout << "Waves Type: " << waveConfig.type << std::endl;
    std::cout << "Wave Height: " << waveConfig.height << std::endl;
    std::cout << "Wave Period: " << waveConfig.period << std::endl;
    std::cout << "Wave Spectrum Type: " << waveConfig.spectrumType << std::endl;
    std::cout << "Wave Directions: ";
    for (double dir : waveConfig.direction) {
        std::cout << dir << " ";
    }
    std::cout << std::endl;
    std::cout << "Wave Elevation File: " << waveConfig.elevationFile << std::endl;

    return waveConfig;
}

std::map<int, BodyConfig> parseBodyConfig(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        throw std::runtime_error("Failed to open file.");
    }

    // Extracting directory from the input file path
    std::filesystem::path inputPath(filename);
    std::string directory = inputPath.parent_path().string();

    std::map<int, BodyConfig> bodies;
    std::string line;
    BodyConfig body;
    int bodyNumber = 0;  // Start at 0, increment when a new body is found

    while (std::getline(inFile, line)) {
        if (line.find("bodyClass('") != std::string::npos) {
            if (bodyNumber > 0) {  // Save previous body if exists
                bodies[bodyNumber] = body;
            }
            bodyNumber++;         // Increment for new body
            body = BodyConfig();  // Reset body

            size_t start       = line.find("'") + 1;
            size_t end         = line.find_last_of("'");
            body.hydroDataFile = line.substr(start, end - start);
            body.hydroDataFile = directory + "/" + body.hydroDataFile;
            std::cout << "HydroDataFile: " << body.hydroDataFile << std::endl;
        } else if (line.find(".geometryFile = '") != std::string::npos) {
            size_t start      = line.find("'") + 1;
            size_t end        = line.find_last_of("'");
            body.geometryFile = line.substr(start, end - start);
            body.geometryFile = directory + "/" + body.geometryFile;
            std::cout << "GeometryFile: " << body.geometryFile << std::endl;
        } else if (line.find(".mass = ") != std::string::npos) {
            size_t start = line.find(".mass = ") + 8;
            body.mass    = std::stod(line.substr(start));
            std::cout << "Mass: " << body.mass << std::endl;
        } else if (line.find(".inertia = [") != std::string::npos) {
            size_t start           = line.find("[") + 1;
            size_t end             = line.find("]");
            std::string inertiaStr = line.substr(start, end - start);
            std::istringstream iss(inertiaStr);
            double val;
            body.inertia.clear();  // Clear any existing values
            while (iss >> val) {
                body.inertia.push_back(val);
                if (iss.peek() == ' ') iss.ignore();
            }
            std::cout << "Inertia: ";
            for (double i : body.inertia)
                std::cout << i << " ";
            std::cout << std::endl;
        } else if (line.find(".position = [") != std::string::npos) {
            size_t start            = line.find("[") + 1;
            size_t end              = line.find("]");
            std::string positionStr = line.substr(start, end - start);

            // Replacing commas with spaces to handle both as delimiters
            std::replace(positionStr.begin(), positionStr.end(), ',', ' ');

            std::istringstream iss(positionStr);
            double val;
            body.position.clear();  // Clear any existing values
            while (iss >> val) {
                body.position.push_back(val);
            }

            std::cout << "Position: ";
            for (double p : body.position)
                std::cout << p << " ";
            std::cout << std::endl;
        }
    }

    // Save the last body
    if (bodyNumber > 0) {
        bodies[bodyNumber] = body;
    }

    return bodies;
}

std::vector<PTOConfig> parsePTOConfig(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        throw std::runtime_error("Failed to open file.");
    }

    std::vector<PTOConfig> configs;
    PTOConfig currentConfig;
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.find("ptoClass('") != std::string::npos) {
            // Save the previous config if it is not empty
            if (!currentConfig.type.empty()) {
                configs.push_back(currentConfig);
                currentConfig = PTOConfig();  // Reset current config
            }

            currentConfig.type = line.substr(line.find("('") + 2, line.find("')") - line.find("('") - 2);
        } else if (line.find(".stiffness = ") != std::string::npos) {
            currentConfig.stiffness = std::stod(line.substr(line.find("=") + 2));
        } else if (line.find(".damping = ") != std::string::npos) {
            currentConfig.damping = std::stod(line.substr(line.find("=") + 2));
        } else if (line.find(".location = [") != std::string::npos) {
            currentConfig.location.clear();  // Ensure the location is cleared before adding new values

            std::string loc = line.substr(line.find("[") + 1, line.find("]") - line.find("[") - 1);
            std::istringstream locStream(loc);
            double val;
            while (locStream >> val) {
                currentConfig.location.push_back(val);
                if (locStream.peek() == ',') {
                    locStream.ignore();
                }
            }
        } else if (line.find(".bodies = [") != std::string::npos) {
            size_t start          = line.find("[") + 1;
            size_t end            = line.find("]");
            std::string bodiesStr = line.substr(start, end - start);

            size_t pos;
            while ((pos = bodiesStr.find("body(")) != std::string::npos) {
                size_t startPos        = pos + 5;  // position after "body("
                size_t endPos          = bodiesStr.find(")", startPos);
                std::string bodyNumStr = bodiesStr.substr(startPos, endPos - startPos);

                int bodyNum = std::stoi(bodyNumStr);
                currentConfig.bodies.push_back(bodyNum);

                bodiesStr = bodiesStr.substr(endPos + 1);  // get the remaining string
            }
        } else if (line.find(".attachments = ") != std::string::npos) {
            size_t start               = line.find("[[") + 2;
            size_t end                 = line.find("]]");
            std::string attachmentsStr = line.substr(start, end - start);

            std::istringstream iss(attachmentsStr);
            std::string attachmentStr;
            while (getline(iss, attachmentStr, ']')) {
                attachmentStr.erase(std::remove(attachmentStr.begin(), attachmentStr.end(), '['), attachmentStr.end());
                attachmentStr.erase(std::remove(attachmentStr.begin(), attachmentStr.end(), ','), attachmentStr.end());

                std::istringstream attachStream(attachmentStr);
                std::vector<double> attachment;
                double val;
                while (attachStream >> val) {
                    attachment.push_back(val);
                }
                currentConfig.attachments.push_back(attachment);
            }
        }
    }
    // Save the last parsed config
    if (!currentConfig.type.empty()) {
        configs.push_back(currentConfig);
    }
    for (const auto& config : configs) {
        std::cout << "PTO Type: " << config.type << std::endl;
        std::cout << "Stiffness: " << config.stiffness << std::endl;
        std::cout << "Damping: " << config.damping << std::endl;
        std::cout << "Location: ";
        for (auto val : config.location) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        std::cout << "Bodies: ";
        for (auto body : config.bodies) {
            std::cout << body << " ";
        }
        std::cout << std::endl;
        std::cout << "Attachments: " << std::endl;
        for (const auto& vec : config.attachments) {
            std::cout << "[";
            for (auto val : vec) {
                std::cout << val << " ";
            }
            std::cout << "]" << std::endl;
        }
        std::cout << std::endl;
    }
    return configs;
}

std::vector<std::shared_ptr<ChBody>> createBodies(const std::map<int, BodyConfig>& bodyConfigs,
                                                  const std::string& inputFilePath,
                                                  ChSystem& system) {
    std::vector<std::shared_ptr<ChBody>> bodies;

    // Extracting directory from the input file path
    std::filesystem::path inputPath(inputFilePath);
    std::string directory = inputPath.parent_path().string();

    for (const auto& [bodyNumber, bodyConfig] : bodyConfigs) {
        std::string fullMeshFilename = bodyConfig.geometryFile;
        auto body                    = chrono_types::make_shared<ChBodyEasyMesh>(fullMeshFilename,  // mesh filename
                                                              0,                 // density
                                                              false,  // do not evaluate mass automatically
                                                              true,  // create visualization asset
                                                              false  // collisions
        );

        std::cout << "Creating Body Number: " << bodyNumber << std::endl;
        std::cout << "Mesh Filename: " << fullMeshFilename << std::endl;

        // Setting the body parameters
        body->SetNameString("body" + std::to_string(bodyNumber));
        std::cout << "Body Name: " << body->GetNameString() << std::endl;

        if (bodyConfig.position.size() == 3) {
            body->SetPos(ChVector<>(bodyConfig.position[0], bodyConfig.position[1], bodyConfig.position[2]));
            std::cout << "Position: (" << bodyConfig.position[0] << ", " << bodyConfig.position[1] << ", "
                      << bodyConfig.position[2] << ")" << std::endl;
        } else {
            std::cerr << "Error: Position vector size is not 3." << std::endl;
        }

        body->SetMass(bodyConfig.mass);
        std::cout << "Mass: " << bodyConfig.mass << std::endl;

        if (bodyConfig.inertia.size() == 3) {
            body->SetInertiaXX(ChVector<>(bodyConfig.inertia[0], bodyConfig.inertia[1], bodyConfig.inertia[2]));
            std::cout << "Inertia: (" << bodyConfig.inertia[0] << ", " << bodyConfig.inertia[1] << ", "
                      << bodyConfig.inertia[2] << ")" << std::endl;
        } else {
            std::cerr << "Error: Inertia vector size is not 3." << std::endl;
        }

        // Adding body to the system
        system.AddBody(body);
        std::cout << "Body added to the system." << std::endl;

        // Adding body to the list of bodies
        bodies.push_back(body);
        std::cout << "Body added to the list of bodies." << std::endl;
        std::cout << std::endl;
    }

    std::cout << "Total bodies created: " << bodies.size() << std::endl;
    std::cout << std::endl;

    return bodies;
}

// Function to create joints or PTOs
void createJointOrPTO(ChSystem& system,
                      const std::string& type,
                      std::shared_ptr<ChBody> body1,
                      std::shared_ptr<ChBody> body2,
                      const ChVector<>& location,
                      const ChVector<>& attachment1,
                      const ChVector<>& attachment2,
                      double stiffness,
                      double damping) {
    if (type == "Translational") {
        // Creating a prismatic joint
        auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
        prismatic->Initialize(body1, body2, false, ChCoordsys<>(attachment1), ChCoordsys<>(attachment2));
        system.AddLink(prismatic);

        // Creating the PTO with translational spring-damper
        auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
        prismatic_pto->Initialize(body1, body2, false, attachment1, attachment2);
        prismatic_pto->SetSpringCoefficient(stiffness);
        prismatic_pto->SetDampingCoefficient(damping);
        system.AddLink(prismatic_pto);
    }
    // else if (type == "Rotational") {
    //     // Code for creating a rotational joint or PTO
    // }
}

std::shared_ptr<WaveBase> setupWaveParameters(const WaveConfig& waveConfig,
                                              double timestep,
                                              double simulationDuration,
                                              double rampTime,
                                              size_t num_bodies) {
    std::shared_ptr<WaveBase> hydro_inputs;

    std::cout << "Setting up wave parameters:" << std::endl;
    std::cout << "Wave Type: " << waveConfig.type << std::endl;
    if (waveConfig.type == "irregular") {
        IrregularWaveParams wave_inputs;
        wave_inputs.num_bodies_          = num_bodies;
        wave_inputs.simulation_dt_       = timestep;
        wave_inputs.simulation_duration_ = simulationDuration;
        wave_inputs.ramp_duration_       = rampTime;
        wave_inputs.wave_height_         = waveConfig.height;
        wave_inputs.wave_period_         = waveConfig.period;

        std::cout << "Simulation Timestep: " << wave_inputs.simulation_dt_ << std::endl;
        std::cout << "Simulation Duration: " << wave_inputs.simulation_duration_ << std::endl;
        std::cout << "Ramp Time: " << wave_inputs.ramp_duration_ << std::endl;
        std::cout << "Wave Height: " << wave_inputs.wave_height_ << std::endl;
        std::cout << "Wave Period: " << wave_inputs.wave_period_ << std::endl;

        try {
            hydro_inputs = std::make_shared<IrregularWaves>(wave_inputs);
        } catch (const std::exception& e) {
            std::cerr << "Caught exception: " << e.what() << '\n';
        } catch (...) {
            std::cerr << "Caught unknown exception.\n";
        }
    } else if (waveConfig.type == "regular") {
        auto regular_wave = std::make_shared<RegularWave>(num_bodies);  // or any other appropriate value
        regular_wave->regular_wave_amplitude_ = waveConfig.height / 2;
        regular_wave->regular_wave_omega_ = 2.0 * CH_C_PI / waveConfig.period;  // calculating omega based on the period

        std::cout << "Wave Amplitude: " << regular_wave->regular_wave_amplitude_ << std::endl;
        std::cout << "Wave Omega: " << regular_wave->regular_wave_omega_ << std::endl;

        hydro_inputs = regular_wave;  // Assigning back to hydro_inputs of type WaveBase
    } else if (waveConfig.type == "still") {
        std::cout << "No wave conditions (Still water)" << std::endl;
        hydro_inputs = std::make_shared<NoWave>(1);
    }

    return hydro_inputs;
}

void initializeOutputFile(std::ofstream& outputFile, const std::vector<std::shared_ptr<ChBody>>& bodies) {
    outputFile.open("./results/output.txt");

    if (!outputFile.is_open()) {
        if (!std::filesystem::exists("./results")) {
            std::filesystem::create_directory("./results");
            outputFile.open("./results/output.txt");
        }
    }

    outputFile << std::left << std::setw(20) << "Time (s)";
    for (int i = 0; i < bodies.size(); ++i) {
        std::string bodyIndex = "Body" + std::to_string(i);
        outputFile << std::setw(16) << (bodyIndex + "_x (m)") << std::setw(16) << (bodyIndex + "_y (m)")
                   << std::setw(16) << (bodyIndex + "_z (m)");
    }
    outputFile << std::endl;
}

void saveDataToFile(std::ofstream& outputFile,
                    const std::vector<double>& time_vector,
                    const std::map<int, std::vector<ChVector<>>>& body_positions) {
    for (int i = 0; i < time_vector.size(); ++i) {
        outputFile << std::left << std::setw(20) << std::setprecision(2) << std::fixed << time_vector[i];

        for (const auto& [bodyIndex, positions] : body_positions) {
            outputFile << std::setw(16) << std::setprecision(4) << std::fixed << positions[i].x() << std::setw(16)
                       << std::setprecision(4) << std::fixed << positions[i].y() << std::setw(16)
                       << std::setprecision(4) << std::fixed << positions[i].z();
        }
        outputFile << std::endl;
    }
    outputFile.close();
}

// usage: ./<demos>.exe [WECSimInputFile.m] [--nogui]

int main(int argc, char* argv[]) {
    GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";

    if (hydroc::SetInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path to .m file>" << std::endl;
        return 1;
    }

    std::string filePath = argv[1];

    // Parsing Simulation Configuration
    SimulationConfig simuConfig = parseSimulationConfig(filePath);

    // Parsing Wave Configuration
    WaveConfig waveConfig = parseWaveConfig(filePath);

    // Parsing Body Information
    std::map<int, BodyConfig> bodyConfigs = parseBodyConfig(filePath);

    // Parsing PTO Configuration
    std::vector<PTOConfig> ptoConfigs = parsePTOConfig(filePath);

    // System/solver settings
    ChSystemNSC system;

    // Defining some standard solver settings
    system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
    system.SetTimestepperType(ChTimestepper::Type::HHT);
    system.SetSolverType(ChSolver::Type::GMRES);
    system.SetSolverMaxIterations(300);
    ChRealtimeStepTimer realtime_timer;

    // Defining simulation duration and timestep
    double simulationDuration = simuConfig.endTime;
    double timestep           = simuConfig.dt;
    system.SetStep(timestep);

    // Create body objects
    std::vector<std::shared_ptr<ChBody>> bodies = createBodies(bodyConfigs, filePath, system);

    for (const auto& ptoConfig : ptoConfigs) {
        // Extracting the information from each PTOConfig
        std::string type              = ptoConfig.type;
        int indexBody1                = ptoConfig.bodies[0] - 1;  // get the index of the first body
        int indexBody2                = ptoConfig.bodies[1] - 1;  // get the index of the second body
        std::shared_ptr<ChBody> body1 = bodies[indexBody1];       // Get the first body using the index
        std::shared_ptr<ChBody> body2 = bodies[indexBody2];       // Get the second body using the index
        ChVector<> location;
        ChVector<> attachment1;
        ChVector<> attachment2;
        if (ptoConfig.location.size() != 3) {
            std::cerr << "Error: ptoConfig.location has " << ptoConfig.location.size() << " elements: ";
            for (const auto& val : ptoConfig.location) {
                std::cerr << val << " ";
            }
            std::cerr << std::endl;
        }
        if (ptoConfig.attachments[0].size() == 3 && ptoConfig.attachments[1].size() == 3) {
            attachment1.Set(ptoConfig.attachments[0][0], ptoConfig.attachments[0][1], ptoConfig.attachments[0][2]);
            attachment2.Set(ptoConfig.attachments[1][0], ptoConfig.attachments[1][1], ptoConfig.attachments[1][2]);
        } else {
            // Handle the error
            std::cerr << "Error: ptoConfig.attachments elements must each have exactly 3 elements." << std::endl;
        }
        double stiffness = ptoConfig.stiffness;
        double damping   = ptoConfig.damping;

        // Calling createJointOrPTO() function with the extracted information
        createJointOrPTO(system, type, body1, body2, location, attachment1, attachment2, stiffness, damping);
    }

    auto hydro_inputs = setupWaveParameters(waveConfig,
                                            simuConfig.dt,
                                            simuConfig.endTime,
                                            simuConfig.rampTime,
                                            bodies.size());

    std::cerr << "hydroDataFile location: " << bodyConfigs[1].hydroDataFile << std::endl;

    TestHydro hydro_forces(bodies, bodyConfigs[1].hydroDataFile);
    hydro_forces.AddWaves(hydro_inputs);

    // for profiling
    auto start = std::chrono::high_resolution_clock::now();

    // for visualization
    bool visualizationOn = true;
    if (argc > 2 && std::string("--nogui").compare(argv[2]) == 0) {
        visualizationOn = false;
    }

    // for output
    std::vector<double> time_vector;
    std::map<int, std::vector<ChVector<>>> body_positions;

    // create an ofstream object
    std::ofstream outputFile;

    // initialize the output file
    initializeOutputFile(outputFile, bodies);

    // main simulation loop
    std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);
    hydroc::gui::UI& ui                  = *pui.get();
    ui.Init(&system, "RM3 - Regular Wave Test");
    ui.SetCamera(0, -50, -10, 0, 0, -10);

    while (system.GetChTime() <= simulationDuration) {
        if (ui.IsRunning(timestep) == false) break;

        if (ui.simulationStarted) {
            system.DoStepDynamics(timestep);
            time_vector.push_back(system.GetChTime());

            for (int i = 0; i < bodies.size(); ++i) {
                body_positions[i].push_back(bodies[i]->GetPos());
            }
        }
    }

    saveDataToFile(outputFile, time_vector, body_positions);

    return 0;
}

//// If no argument is given user can set HYDROCHRONO_DATA_DIR
//// environment variable to give the data_directory.
////
// int main(int argc, char* argv[]) {
//    GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";
//
//    if (hydroc::SetInitialEnvironment(argc, argv) != 0) {
//        return 1;
//    }
//
//    // Check if --nogui option is set as 2nd argument
//    bool visualizationOn = true;
//    if (argc > 2 && std::string("--nogui").compare(argv[2]) == 0) {
//        visualizationOn = false;
//    }
//
//    // Get model file names
//    std::filesystem::path DATADIR(hydroc::getDataDir());
//
//    auto body1_meshfame = (DATADIR / "rm3" / "geometry" / "float_cog.obj").lexically_normal().generic_string();
//    auto body2_meshfame = (DATADIR / "rm3" / "geometry" / "plate_cog.obj").lexically_normal().generic_string();
//    auto h5fname        = (DATADIR / "rm3" / "hydroData" / "rm3.h5").lexically_normal().generic_string();
//
//    // system/solver settings
//    ChSystemNSC system;
//
//    system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
//    double timestep = 0.01;
//    system.SetTimestepperType(ChTimestepper::Type::HHT);
//    system.SetSolverType(ChSolver::Type::GMRES);
//    system.SetSolverMaxIterations(300);  // the higher, the easier to keep the constraints satisfied.
//    system.SetStep(timestep);
//    ChRealtimeStepTimer realtime_timer;
//    double simulationDuration = 40.0;
//
//    // Create user interface
//    std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);
//
//    hydroc::gui::UI& ui = *pui.get();
//
//    // some io/viz options
//    bool profilingOn = true;
//    bool saveDataOn  = true;
//    std::vector<double> time_vector;
//    std::vector<double> float_heave_position;
//    std::vector<double> float_drift_position;
//    std::vector<double> plate_heave_position;
//
//    // set up body from a mesh
//    std::cout << "Attempting to open mesh file: " << body1_meshfame << std::endl;
//    std::shared_ptr<ChBody> float_body1 = chrono_types::make_shared<ChBodyEasyMesh>(  //
//        body1_meshfame,
//        0,      // density
//        false,  // do not evaluate mass automatically
//        true,   // create visualization asset
//        false   // collisions
//    );
//
//    // define the float's initial conditions
//    system.Add(float_body1);
//    float_body1->SetNameString("body1");
//    float_body1->SetPos(ChVector<>(0, 0, -0.72));
//    float_body1->SetMass(725834);
//    float_body1->SetInertiaXX(ChVector<>(20907301.0, 21306090.66, 37085481.11));
//    // float_body1->SetCollide(false);
//
//    // Create a visualization material
//    auto red = chrono_types::make_shared<ChVisualMaterial>();
//    red->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.1f));
//    float_body1->GetVisualShape(0)->SetMaterial(0, red);
//
//    std::cout << "Attempting to open mesh file: " << body2_meshfame << std::endl;
//    std::shared_ptr<ChBody> plate_body2 = chrono_types::make_shared<ChBodyEasyMesh>(  //
//        body2_meshfame,
//        0,      // density
//        false,  // do not evaluate mass automatically
//        true,   // create visualization asset
//        false   // collisions
//    );
//
//    // Create a visualization material
//    auto blue = chrono_types::make_shared<ChVisualMaterial>();
//    blue->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.6f));
//    plate_body2->GetVisualShape(0)->SetMaterial(0, blue);
//
//    // define the plate's initial conditions
//    system.Add(plate_body2);
//    plate_body2->SetNameString("body2");
//    plate_body2->SetPos(ChVector<>(0, 0, (-21.29)));
//    plate_body2->SetMass(886691);
//    plate_body2->SetInertiaXX(ChVector<>(94419614.57, 94407091.24, 28542224.82));
//    // plate_body2->SetCollide(false);
//
//    // add prismatic joint between the two bodies
//    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
//    prismatic->Initialize(float_body1, plate_body2, false, ChCoordsys<>(ChVector<>(0, 0, -0.72)),
//                          ChCoordsys<>(ChVector<>(0, 0, -21.29)));
//    system.AddLink(prismatic);
//
//    auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
//    prismatic_pto->Initialize(float_body1, plate_body2, false, ChVector<>(0, 0, -0.72), ChVector<>(0, 0, -21.29));
//    prismatic_pto->SetDampingCoefficient(0.0);
//    system.AddLink(prismatic_pto);
//
//    // define wave parameters
//    auto my_hydro_inputs                    = std::make_shared<RegularWave>();
//    my_hydro_inputs->regular_wave_amplitude_ = 1.0;
//    my_hydro_inputs->regular_wave_omega_     = 2.10;
//
//    // attach hydrodynamic forces to body
//    std::vector<std::shared_ptr<ChBody>> bodies;
//    bodies.push_back(float_body1);
//    bodies.push_back(plate_body2);
//    TestHydro hydro_forces(bodies, h5fname);
//    hydro_forces.AddWaves(my_hydro_inputs);
//
//    // for profiling
//    auto start = std::chrono::high_resolution_clock::now();
//
//    // main simulation loop
//    ui.Init(&system, "RM3 - Regular Wave Test");
//    ui.SetCamera(0, -50, -10, 0, 0, -10);
//
//    while (system.GetChTime() <= simulationDuration) {
//        if (ui.IsRunning(timestep) == false) break;
//
//        if (ui.simulationStarted) {
//            system.DoStepDynamics(timestep);
//
//            // append data to output vector
//            time_vector.push_back(system.GetChTime());
//            float_heave_position.push_back(float_body1->GetPos().z());
//            float_drift_position.push_back(float_body1->GetPos().x());
//            plate_heave_position.push_back(plate_body2->GetPos().z());
//        }
//    }
//
//    // for profiling
//    auto end          = std::chrono::high_resolution_clock::now();
//    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//
//    if (profilingOn) {
//        std::ofstream profilingFile;
//        profilingFile.open("./results/rm3_reg_waves_duration.txt");
//        if (!profilingFile.is_open()) {
//            if (!std::filesystem::exists("./results")) {
//                std::cout << "Path " << std::filesystem::absolute("./results")
//                          << " does not exist, creating it now..." << std::endl;
//                std::filesystem::create_directory("./results");
//                profilingFile.open("./results/rm3_reg_waves_duration.txt");
//                if (!profilingFile.is_open()) {
//                    std::cout << "Still cannot open file, ending program" << std::endl;
//                    return 0;
//                }
//            }
//        }
//        profilingFile << duration << " ms\n";
//        profilingFile.close();
//    }
//
//    if (saveDataOn) {
//        std::ofstream outputFile;
//        outputFile.open("./results/rm3_reg_waves.txt");
//        if (!outputFile.is_open()) {
//            if (!std::filesystem::exists("./results")) {
//                std::cout << "Path " << std::filesystem::absolute("./results")
//                          << " does not exist, creating it now..." << std::endl;
//                std::filesystem::create_directory("./results");
//                outputFile.open("./results/rm3_decay.txt");
//                if (!outputFile.is_open()) {
//                    std::cout << "Still cannot open file, ending program" << std::endl;
//                    return 0;
//                }
//            }
//        }
//        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(16) << "Float Heave (m)"
//                   << std::right << std::setw(16) << "Plate Heave (m)" << std::right << std::setw(16)
//                   << "Float Drift (x) (m)" << std::endl;
//        for (int i = 0; i < time_vector.size(); ++i)
//            outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
//                       << std::right << std::setw(16) << std::setprecision(4) << std::fixed << float_heave_position[i]
//                       << std::right << std::setw(16) << std::setprecision(4) << std::fixed << plate_heave_position[i]
//                       << std::right << std::setw(16) << std::setprecision(4) << std::fixed << float_drift_position[i]
//                       << std::endl;
//        outputFile.close();
//    }
//    return 0;
//}