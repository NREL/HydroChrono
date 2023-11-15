// Include necessary headers
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

// Use namespaces for cleaner code and to avoid name conflicts
using namespace chrono;
using namespace chrono::geometry;

// Define a struct for simulation configurations
struct SimulationConfig {
    std::string explorer;
    std::string modelName;
    double startTime;
    double rampTime;
    double endTime;
    double dt;  // time step
};

// Define a struct for wave configurations
struct WaveConfig {
    std::string type;               // type of wave
    double height;                  // height of the wave
    double period;                  // period of the wave
    std::string spectrumType;       // type of spectrum used
    std::vector<double> direction;  // direction of the wave
    std::string elevationFile;      // file containing elevation data
};

// Define a struct for body configurations
struct BodyConfig {
    std::string hydroDataFile;  // file containing hydrodynamic data
    std::string type;              // Type of the body
    double mass;                   // Mass of the body
    std::vector<double> inertia;   // Inertia of the body
    std::vector<double> position;  // Position of the body
    // For mesh bodies
    std::string geometryFile;  // File path for the geometry
    // For sphere bodies
    double radius;   // Radius of the sphere
    // ... other properties as needed
};

// Define a struct for PTO configurations
struct PTOConfig {
    std::string type;                              // Type of PTO
    std::vector<int> bodies;                       // Bodies connected to the PTO
    std::vector<std::vector<double>> attachments;  // Attachment points
    double stiffness;                              // Stiffness of the PTO (used in spring-damper, prismatic, etc.)
    double damping;                                // Damping factor (used in spring-damper, prismatic, etc.)

    // Optional fields for specific joint types
    // For spring-damper
    double rest_length;   // Rest length of the spring
    double spring_radius  = 0.5;  // Radius of the spring (for visualization)
    int spring_resolution = 1280;    // Resolution of the spring shape
    int spring_turns      = 16;    // Number of turns in the spring shape

    // Additional fields can be added here for other joint types
};

// Function to parse the Simulation Configuration from a given file
SimulationConfig parseSimulationConfig(const std::string& filename) {
    // Open the file for reading
    std::ifstream inFile(filename);
    if (!inFile) {
        // Throw an error if the file can't be opened
        throw std::runtime_error("Failed to open file: " + filename);
    }

    SimulationConfig config;
    std::string line;
    while (std::getline(inFile, line)) {  // Read each line from the file
        std::istringstream iss(line);
        std::string token;

        // Parse and assign values to the corresponding attributes of the SimulationConfig struct
        if (line.find("simu.modelName = '") != std::string::npos) {
            config.modelName = line.substr(line.find("'") + 1, line.rfind("'") - line.find("'") - 1);
        } else if (line.find("simu.explorer = '") != std::string::npos) {
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

    // Output the parsed values to the console for verification
    //std::cout << "Parsed Simulation Configuration:" << std::endl;
    //std::cout << "Model Name: " << config.modelName << std::endl;
    //std::cout << "Explorer: " << config.explorer << std::endl;
    //std::cout << "Start Time: " << config.startTime << std::endl;
    //std::cout << "Ramp Time: " << config.rampTime << std::endl;
    //std::cout << "End Time: " << config.endTime << std::endl;
    //std::cout << "Time Step (dt): " << config.dt << std::endl;

    // Return the filled SimulationConfig struct
    return config;
}

// Function to parse the Wave Configuration from a given file
WaveConfig parseWaveConfig(const std::string& filePath) {
    std::ifstream file(filePath);
    std::string line;
    WaveConfig waveConfig;

    // Check if the file is successfully opened
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return waveConfig;
    }

    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string word;
        iss >> word;

        // Parse and assign values to the corresponding attributes of the WaveConfig struct
        if (word == "waves") {
            iss >> word;  // Skip '='
            std::string temp;
            iss >> temp;  // Getting the full string, e.g., waveClass('regular');
            // Extracting the word within single quotes
            size_t start    = temp.find("'") + 1;
            size_t end      = temp.find_last_of("'");
            waveConfig.type = temp.substr(start, end - start);
        } else if (word == "waves.height") {
            iss >> word;  // Skip '='
            iss >> waveConfig.height;
        } else if (word == "waves.period") {
            iss >> word;  // Skip '='
            iss >> waveConfig.period;
        } else if (word == "waves.spectrumType") {
            iss >> word;  // Skip '='
            iss >> waveConfig.spectrumType;
        } else if (word == "waves.direction") {
            iss >> word;  // Skip '='
            double direction;
            while (iss >> direction) {
                waveConfig.direction.push_back(direction);
            }
        } else if (word == "waves.elevationFile") {
            iss >> word;  // Skip '='
            iss >> waveConfig.elevationFile;
        }
    }

    file.close();

    //// Output the parsed values to the console for verification
    std::cout << "Parsed Wave Configuration:" << std::endl;
    std::cout << "Waves Type: " << waveConfig.type << std::endl;
    //std::cout << "Wave Height: " << waveConfig.height << std::endl;
    //std::cout << "Wave Period: " << waveConfig.period << std::endl;
    //std::cout << "Wave Spectrum Type: " << waveConfig.spectrumType << std::endl;
    //std::cout << "Wave Directions: ";
    //for (double dir : waveConfig.direction) {
    //    std::cout << dir << " ";
    //}
    //std::cout << std::endl;
    //std::cout << "Wave Elevation File: " << waveConfig.elevationFile << std::endl;

    // Return the filled WaveConfig struct
    return waveConfig;
}

// Helper function to extract strings within single quotes and assign them
void extractAndAssignString(const std::string& line,
                            std::string& attribute,
                            const std::string& directory,
                            bool appendDirectory = true) {
    size_t start = line.find("'") + 1;
    size_t end   = line.find_last_of("'");
    attribute    = line.substr(start, end - start);
    if (appendDirectory) {
        attribute = directory + "/" + attribute;
    }
}
// Helper function to extract and assign double values
double extractAndAssignDouble(const std::string& line, const std::string& pattern) {
    size_t start = line.find(pattern) + pattern.length();
    return std::stod(line.substr(start));
}

// Helper function to extract and assign vectors
void extractAndAssignVector(const std::string& line, std::vector<double>& attribute) {
    size_t start        = line.find("[") + 1;
    size_t end          = line.find("]");
    std::string numbers = line.substr(start, end - start);

    // Replacing commas with spaces to handle both as delimiters
    std::replace(numbers.begin(), numbers.end(), ',', ' ');

    std::istringstream iss(numbers);
    double val;
    attribute.clear();  // Clear any existing values
    while (iss >> val) {
        attribute.push_back(val);
    }
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

            // Extract and assign hydroDataFile
            extractAndAssignString(line, body.hydroDataFile, directory);
            //std::cout << "HydroDataFile: " << body.hydroDataFile << std::endl;

        } else if (line.find(".type = '") != std::string::npos) {
            extractAndAssignString(line, body.type, directory, false);
        } else if (line.find(".radius = ") != std::string::npos) {
            body.radius = extractAndAssignDouble(line, ".radius = ");
        } else if (line.find(".geometryFile = '") != std::string::npos) {
            // Extract and assign geometryFile
            extractAndAssignString(line, body.geometryFile, directory);
            //std::cout << "GeometryFile: " << body.geometryFile << std::endl;
        } else if (line.find(".mass = ") != std::string::npos) {
            body.mass = extractAndAssignDouble(line, ".mass = ");
            //std::cout << "Mass: " << body.mass << std::endl;

        } else if (line.find(".inertia = [") != std::string::npos) {
            extractAndAssignVector(line, body.inertia);
            //std::cout << "Inertia: ";
            //for (double i : body.inertia)
            //    std::cout << i << " ";
            //std::cout << std::endl;

        } else if (line.find(".position = [") != std::string::npos) {
            extractAndAssignVector(line, body.position);
            //std::cout << "Position: ";
            //for (double p : body.position)
            //    std::cout << p << " ";
            //std::cout << std::endl;
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
        } else if (line.find(".restLength = ") != std::string::npos) {
            try {
                currentConfig.rest_length = std::stod(line.substr(line.find("=") + 2));
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument: " << ia.what() << " in parsing restLength" << std::endl;
            }
        }
    }
    // Save the last parsed config
    if (!currentConfig.type.empty()) {
        configs.push_back(currentConfig);
    }
    if (configs.empty()) {
        std::cout << "No PTO configurations found in the file." << std::endl;
    } else {
        for (const auto& config : configs) {
            std::cout << "PTO Type: " << config.type << std::endl;
            //    std::cout << "Stiffness: " << config.stiffness << std::endl;
            //    std::cout << "Damping: " << config.damping << std::endl;
            //    std::cout << "Location: ";
            //    for (auto val : config.location) {
            //        std::cout << val << " ";
            //    }
            //    std::cout << std::endl;
            //    std::cout << "Bodies: ";
            //    for (auto body : config.bodies) {
            //        std::cout << body << " ";
            //    }
            //    std::cout << std::endl;
            //    std::cout << "Attachments: " << std::endl;
            //    for (const auto& vec : config.attachments) {
            //        std::cout << "[";
            //        for (auto val : vec) {
            //            std::cout << val << " ";
            //        }
            //        std::cout << "]" << std::endl;
            //    }
            //    std::cout << std::endl;
        }
    }
    return configs;
}

// Set up bodies
void setBodyPosition(std::shared_ptr<ChBodyEasyMesh>& body, const std::vector<double>& position) {
    if (position.size() == 3) {
        body->SetPos(ChVector<>(position[0], position[1], position[2]));
    } else {
        std::cerr << "Error: Position vector size is not 3." << std::endl;
    }
}

void setBodyInertia(std::shared_ptr<ChBodyEasyMesh>& body, const std::vector<double>& inertia) {
    if (inertia.size() == 3) {
        body->SetInertiaXX(ChVector<>(inertia[0], inertia[1], inertia[2]));
    } else {
        std::cerr << "Error: Inertia vector size is not 3." << std::endl;
    }
}

std::shared_ptr<ChBody> createSphereBody(int bodyNumber, const BodyConfig& bodyConfig, ChSystem& system) {
    auto irm = chrono_types::make_shared<ChBodyEasySphere>(bodyConfig.radius, 1000.0, true, false);
    if (!irm) {
        return nullptr;
    }

    irm->SetMass(bodyConfig.mass);
    irm->SetInertiaXX(ChVector<>(bodyConfig.inertia[0], bodyConfig.inertia[1], bodyConfig.inertia[2]));
    irm->SetPos(ChVector<>(bodyConfig.position[0], bodyConfig.position[1], bodyConfig.position[2]));
    irm->SetNameString("body" + std::to_string(bodyNumber));

    system.AddBody(irm);
    return irm;
}

std::shared_ptr<ChBody> createSingleBody(int bodyNumber, const BodyConfig& bodyConfig, ChSystem& system) {
    // Debug print to check the geometry file path
    std::cout << "Creating body " << bodyNumber << std::endl;
    std::cout << "Geometry file: " << bodyConfig.geometryFile << std::endl;

    auto body = chrono_types::make_shared<ChBodyEasyMesh>(bodyConfig.geometryFile, 0, false, true, false);

    if (!body) {
        std::cerr << "Failed to create body " << bodyNumber << std::endl;
        return nullptr;
    }

    body->SetNameString("body" + std::to_string(bodyNumber));
    setBodyPosition(body, bodyConfig.position);
    body->SetMass(bodyConfig.mass);
    setBodyInertia(body, bodyConfig.inertia);

    system.AddBody(body);
    return body;
}

std::pair<std::vector<std::shared_ptr<ChBody>>, std::vector<std::shared_ptr<ChBody>>>
createBodies(const std::map<int, BodyConfig>& bodyConfigs, const std::string& inputFilePath, ChSystem& system) {
    std::vector<std::shared_ptr<ChBody>> allBodies;
    std::vector<std::shared_ptr<ChBody>> hydrodynamicBodies;

    for (const auto& [bodyNumber, bodyConfig] : bodyConfigs) {
        std::shared_ptr<ChBody> body;

        std::cout << "Creating body " << bodyNumber << std::endl;
        std::cout << "Body type: " << bodyConfig.type << std::endl;

        if (bodyConfig.type == "nonhydro-sphere") {
            std::cout << "Creating nonhydro-sphere body." << std::endl;
            body = createSphereBody(bodyNumber, bodyConfig, system);
        } else {
            std::cout << "Creating hydrodynamic body with geometry file: " << bodyConfig.geometryFile << std::endl;
            body = createSingleBody(bodyNumber, bodyConfig, system);
            if (body) {
                hydrodynamicBodies.push_back(body);
            }
        }

        if (body) {
            allBodies.push_back(body);
            std::cout << "Body " << bodyNumber << " added successfully." << std::endl;
        } else {
            std::cerr << "Error creating body " << bodyNumber << "." << std::endl;
        }
    }

    std::cout << "Total bodies created: " << allBodies.size() << std::endl;
    std::cout << "Total hydrodynamic bodies created: " << hydrodynamicBodies.size() << std::endl;

    return std::make_pair(allBodies, hydrodynamicBodies);
}

// Functions to create joints or PTOs
void CreateTranslationalPTO(ChSystem& system,
                              const std::vector<std::shared_ptr<ChBody>>& bodies,
                              const PTOConfig& pto_config) {
    auto body1 = bodies[pto_config.bodies[0] - 1];
    auto body2 = bodies[pto_config.bodies[1] - 1];
    ChVector<> attachment1(pto_config.attachments[0][0], pto_config.attachments[0][1], pto_config.attachments[0][2]);
    ChVector<> attachment2(pto_config.attachments[1][0], pto_config.attachments[1][1], pto_config.attachments[1][2]);

    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(body1, body2, false, ChCoordsys<>(attachment1), ChCoordsys<>(attachment2));
    system.AddLink(prismatic);

    auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
    prismatic_pto->Initialize(body1, body2, false, attachment1, attachment2);
    prismatic_pto->SetSpringCoefficient(pto_config.stiffness);
    prismatic_pto->SetDampingCoefficient(pto_config.damping);
    system.AddLink(prismatic_pto);
}

// Function to create a simple spring-damper
void CreateLinSpringDamper(ChSystem& system,
                        const std::vector<std::shared_ptr<ChBody>>& bodies,
                        const PTOConfig& pto_config,
                        std::vector<std::shared_ptr<ChLinkTSDA>>& ptos) {
    auto body1 = bodies[pto_config.bodies[0] - 1];
    auto body2 = bodies[pto_config.bodies[1] - 1];
    ChVector<> attachment1(pto_config.attachments[0][0], pto_config.attachments[0][1], pto_config.attachments[0][2]);
    ChVector<> attachment2(pto_config.attachments[1][0], pto_config.attachments[1][1], pto_config.attachments[1][2]);

    auto spring_damper = chrono_types::make_shared<ChLinkTSDA>();
    spring_damper->Initialize(body1, body2, true, attachment1, attachment2);
    spring_damper->SetRestLength(pto_config.rest_length);  // Assuming restLength is defined in PTOConfig
    spring_damper->SetSpringCoefficient(pto_config.stiffness);
    spring_damper->SetDampingCoefficient(pto_config.damping);
    spring_damper->AddVisualShape(chrono_types::make_shared<ChSpringShape>(
        pto_config.spring_radius, pto_config.spring_resolution, pto_config.spring_turns));
    system.AddLink(spring_damper);
    ptos.push_back(spring_damper);
}

void CreateJointOrPTO(ChSystem& system,
                      const std::vector<std::shared_ptr<ChBody>>& bodies,
                      const PTOConfig& pto_config,
                      std::vector<std::shared_ptr<ChLinkTSDA>>& ptos) {  // Added ptos vector
    if (pto_config.type == "TranslationalPTO") {
        CreateTranslationalPTO(system, bodies, pto_config);
    } else if (pto_config.type == "LinSpringDamper") {
        CreateLinSpringDamper(system, bodies, pto_config, ptos);  // Pass ptos vector to the function
    }
    // Add more else-if clauses for other joint types
}

std::shared_ptr<WaveBase> setupWaveParameters(const WaveConfig& waveConfig,
                                              double timestep,
                                              double simulationDuration,
                                              double rampTime,
                                              size_t num_bodies) {
    std::shared_ptr<WaveBase> hydro_inputs;

    std::cout << "\n---- Setting up Wave Parameters ----\n";
    std::cout << std::left << std::setw(20) << "Wave Type: " << waveConfig.type << "\n";

    if (waveConfig.type == "irregular") {
        IrregularWaveParams wave_inputs;
        wave_inputs.num_bodies_          = num_bodies;
        wave_inputs.simulation_dt_       = timestep;
        wave_inputs.simulation_duration_ = simulationDuration;
        wave_inputs.ramp_duration_       = rampTime;
        wave_inputs.wave_height_         = waveConfig.height;
        wave_inputs.wave_period_         = waveConfig.period;

        std::cout << std::left << std::setw(20) << "Simulation Timestep:" << wave_inputs.simulation_dt_ << "\n";
        std::cout << std::left << std::setw(20) << "Simulation Duration:" << wave_inputs.simulation_duration_ << "\n";
        std::cout << std::left << std::setw(20) << "Ramp Time:" << wave_inputs.ramp_duration_ << "\n";
        std::cout << std::left << std::setw(20) << "Wave Height:" << wave_inputs.wave_height_ << "\n";
        std::cout << std::left << std::setw(20) << "Wave Period:" << wave_inputs.wave_period_ << "\n";

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

        std::cout << std::left << std::setw(20) << "Wave Amplitude:" << regular_wave->regular_wave_amplitude_ << "\n";
        std::cout << std::left << std::setw(20) << "Wave Omega:" << regular_wave->regular_wave_omega_ << "\n";

        hydro_inputs = regular_wave;  // Assigning back to hydro_inputs of type WaveBase
    } else if (waveConfig.type == "still" || waveConfig.type == "noWaveCIC") {
        hydro_inputs = std::make_shared<NoWave>(num_bodies);
    }
    std::cout << "------------------------------------\n" << std::endl;
    return hydro_inputs;
}

bool openOutputFile(std::ofstream& outputFile, const std::string& filePath) {
    std::filesystem::path outputPath(filePath);

    // Try to open the file directly
    outputFile.open(outputPath);

    // If the file could not be opened, try creating the directory and then opening the file
    if (!outputFile.is_open()) {
        std::filesystem::path dirPath = outputPath.parent_path();
        if (!std::filesystem::exists(dirPath)) {
            std::filesystem::create_directory(dirPath);
        }
        outputFile.open(outputPath);
    }

    return outputFile.is_open();
}

void initializeBodyOutputFile(std::ofstream& bodyOutputFile,
                              const std::vector<std::shared_ptr<ChBody>>& bodies,
                              const std::string& outputDirectory) {
    std::filesystem::path outputDirPath = std::filesystem::path(outputDirectory);

    // Create the directory if it does not exist
    if (!std::filesystem::exists(outputDirPath)) {
        std::filesystem::create_directories(outputDirPath);
    }

    std::filesystem::path outputFile = outputDirPath / "body_output.txt";

    if (openOutputFile(bodyOutputFile, outputFile.string())) {
        bodyOutputFile << std::left << std::setw(20) << "Time (s)";
        for (size_t i = 0; i < bodies.size(); ++i) {
            std::string bodyIndex = "Body" + std::to_string(i);
            bodyOutputFile << std::setw(16) << (bodyIndex + "_x (m)") << std::setw(16) << (bodyIndex + "_y (m)")
                           << std::setw(16) << (bodyIndex + "_z (m)");
        }
        bodyOutputFile << std::endl;
    } else {
        std::cerr << "Failed to open the output file at: " << outputFile << std::endl;
    }
}

void initializePTOOutputFile(std::ofstream& ptoOutputFile, size_t ptoCount, const std::string& outputDirectory) {
    std::filesystem::path outputDirPath = std::filesystem::path(outputDirectory);

    // Create the directory if it does not exist
    if (!std::filesystem::exists(outputDirPath)) {
        std::filesystem::create_directories(outputDirPath);
    }

    std::filesystem::path outputFile = outputDirPath / "pto_output.txt";

    if (openOutputFile(ptoOutputFile, outputFile.string())) {
        ptoOutputFile << std::left << std::setw(20) << "Time (s)";
        for (size_t i = 0; i < ptoCount; ++i) {
            ptoOutputFile << std::setw(20) << ("PTO" + std::to_string(i) + "_power (W)");  // Adjusted column width
        }
        ptoOutputFile << std::endl;
    } else {
        std::cerr << "Failed to open the PTO output file at: " << outputFile << std::endl;
    }
}


void collectPTOData(const std::vector<std::shared_ptr<ChLinkTSDA>>& ptos,
                    std::vector<double>& ptoVelocities,
                    std::vector<double>& ptoPowers) {
    ptoVelocities.clear();
    ptoPowers.clear();

    for (const auto& pto : ptos) {
        double ptoVelocity   = pto->GetVelocity();
        double ptoDamping    = pto->GetDampingCoefficient();
        double powerAbsorbed = ptoVelocity * ptoVelocity * ptoDamping;  // Assuming lower_pto_damping is accessible

        ptoVelocities.push_back(ptoVelocity);
        ptoPowers.push_back(powerAbsorbed);
    }
}

void saveBodyDataToFile(std::ofstream& outputFile,
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

void savePTODataToFile(std::ofstream& ptoOutputFile,
                       const std::vector<double>& time_vector,
                       //const std::vector<std::vector<double>>& ptoVelocities,
                       const std::vector<std::vector<double>>& ptoPowers) {
    for (size_t i = 0; i < time_vector.size(); ++i) {
        ptoOutputFile << std::left << std::setw(20) << std::setprecision(2) << std::fixed << time_vector[i];

        // Assuming each inner vector in ptoVelocities and ptoPowers corresponds to a specific PTO
        for (size_t j = 0; j < ptoPowers.size(); ++j) {
            //if (i < ptoVelocities[j].size()) {
            //    ptoOutputFile << std::setw(16) << std::setprecision(4) << std::fixed << ptoVelocities[j][i];
            //}
            if (i < ptoPowers[j].size()) {
                ptoOutputFile << std::setw(20) << std::setprecision(4) << std::fixed << ptoPowers[j][i];
            }
        }
        ptoOutputFile << std::endl;
    }
    ptoOutputFile.close();
}

// usage: ./<demos>.exe [WECSimInputFile.m] [resultsDirectory] [--nogui]

int main(int argc, char* argv[]) {
    std::cout << R"(
 /\/\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\
         _   _           _            ____ _                                    
        | | | |_   _  __| |_ __ ___  / ___| |__  _ __ ___  _ __   ___                  
        | |_| | | | |/ _` | '__/ _ \| |   | '_ \| '__/ _ \| '_ \ / _ \                 
        |  _  | |_| | (_| | | | (_) | |___| | | | | | (_) | | | | (_) |                
        |_| |_|\__, |\__,_|_|  \___/ \____|_| |_|_|  \___/|_| |_|\___/                 
               |___/                                                            
 /\/\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\
)" << std::endl;

    GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";

    if (hydroc::SetInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <path to .m file> <output directory>" << std::endl;
        return 1;
    }

    std::string filePath        = argv[1];
    std::string outputDirectory = argv[2];

    std::cout << "==== Parsing input file... ====\n" << std::endl;
    
    // Parsing simulation configuration
    SimulationConfig simuConfig = parseSimulationConfig(filePath);

    std::cout << "Model name: " << simuConfig.modelName.c_str() << std::endl;

    // Parsing wave configuration
    std::cout << "Parsing wave config..." << std::endl;
    WaveConfig waveConfig = parseWaveConfig(filePath);

    // Parsing body information
    std::cout << "Parsing body config..." << std::endl;
    std::map<int, BodyConfig> bodyConfigs = parseBodyConfig(filePath);

    // Parsing PTO configuration
    std::cout << "Parsing pto config..." << std::endl;
    std::vector<PTOConfig> ptoConfigs = parsePTOConfig(filePath);

    std::cout << "\n==== Setting up system...  ====" << std::endl;
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
    std::cout << "Creating body objects..." << std::endl;
    auto [bodies, hydrodynamicBodies] = createBodies(bodyConfigs, filePath, system);

    std::cout << "Creating joint and/or PTO objects..." << std::endl;
    std::vector<std::shared_ptr<ChLinkTSDA>> ptos;
    for (const auto& ptoConfig : ptoConfigs) {
        // Validate the attachments and locations
        if (ptoConfig.attachments.size() != 2 || ptoConfig.attachments[0].size() != 3 ||
            ptoConfig.attachments[1].size() != 3) {
            std::cerr << "Error: ptoConfig.attachments elements must each have exactly 3 elements." << std::endl;
            continue;  // Skip this iteration if the data is invalid
        }

        if (ptoConfig.bodies.size() != 2) {
            std::cerr << "Error: ptoConfig.bodies must contain exactly two body indices." << std::endl;
            continue;  // Skip this iteration if the data is invalid
        }

        // Calling createJointOrPTO() function with the PTOConfig object
        CreateJointOrPTO(system, bodies, ptoConfig, ptos);
    }

    std::cout << "Set wave parameters..." << std::endl;
    auto hydro_inputs = setupWaveParameters(waveConfig, simuConfig.dt, simuConfig.endTime, simuConfig.rampTime,
                                            hydrodynamicBodies.size());

    std::cerr << "hydroDataFile location: " << bodyConfigs[1].hydroDataFile << std::endl;

    std::cout << "Test hydro..." << std::endl;
    TestHydro hydro_forces(hydrodynamicBodies, bodyConfigs[1].hydroDataFile);
    hydro_forces.AddWaves(hydro_inputs);

    // For profiling
    auto start = std::chrono::high_resolution_clock::now();

    // For output
    std::vector<double> time;
    std::map<int, std::vector<ChVector<>>> body_positions;
    std::vector<std::vector<double>> ptoVelocities(ptos.size(), std::vector<double>());
    std::vector<std::vector<double>> ptoPowers(ptos.size(), std::vector<double>());

    // Create an ofstream object
    std::ofstream bodyOutputFile;
    std::ofstream ptoOutputFile;

    // Initialize the output file
    std::cout << "Initialize output files..." << std::endl;
    initializeBodyOutputFile(bodyOutputFile, bodies, outputDirectory);
    initializePTOOutputFile(ptoOutputFile, ptoConfigs.size(), outputDirectory);

    // For visualization
    std::cout << std::endl;
    std::cout << "\n---- Initializing Chrono visualization... ----" << std::endl;
    std::cout << std::endl;

    bool visualizationOn = (simuConfig.explorer == "on");
    // Check for command line arguments to override
    if (argc > 3) {
        if (std::string("--nogui").compare(argv[3]) == 0) {
            visualizationOn = false;
        } else if (std::string("--gui").compare(argv[3]) == 0) {
            visualizationOn = true;
        }
    }
    std::cout << "Visualization: " << (visualizationOn ? "Enabled" : "Disabled") << std::endl;

    std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);
    hydroc::gui::UI& ui                  = *pui.get();
    ui.Init(&system, simuConfig.modelName.c_str());
    ui.SetCamera(0, -50, -10, 0, 0, -10);

    // Main simulation loop
    std::cout << std::endl;
    std::cout << "==== Running simulation... ====" << std::endl;
    std::cout << std::endl;

    while (system.GetChTime() <= simulationDuration) {
        if (!ui.IsRunning(timestep)) break;

        if (ui.simulationStarted) {
            system.DoStepDynamics(timestep);
            time.push_back(system.GetChTime());

            // Collect body position data
            for (size_t i = 0; i < bodies.size(); ++i) {
                body_positions[i].push_back(bodies[i]->GetPos());
            }

            // Collect PTO data
            for (size_t i = 0; i < ptos.size(); ++i) {
                double velocity = ptos[i]->GetVelocity();
                double pto_damping_coefficient = ptos[i]->GetDampingCoefficient();
                double power = velocity * velocity * pto_damping_coefficient;  // Assuming pto_damping is accessible
                ptoVelocities[i].push_back(velocity);
                ptoPowers[i].push_back(power);
            }
        }
    }

    saveBodyDataToFile(bodyOutputFile, time, body_positions);
    savePTODataToFile(ptoOutputFile, time, ptoPowers);

    return 0;
}

//// If no argument is given user can set HYDROCHRONO_DATA_DIR
//// environment variable to give the data_directory.
//
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
//    // Create a visualization material
//    auto red = chrono_types::make_shared<ChVisualMaterial>();
//    red->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.1f));
//    float_body1->GetVisualShape(0)->SetMaterial(0, red);
//
//    // Create a visualization material
//    auto blue = chrono_types::make_shared<ChVisualMaterial>();
//    blue->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.6f));
//    plate_body2->GetVisualShape(0)->SetMaterial(0, blue);
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