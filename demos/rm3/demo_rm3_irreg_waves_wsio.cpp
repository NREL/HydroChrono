// Include necessary headers
#include <hydroc/gui/guihelper.h>
#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/core/ChRealtimeStep.h>
//#include <chrono/physics/ChLinkLock.h>
#include <chrono/physics/ChLinkMate.h>  // fixed body uses link

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
    std::string hydroDataFile;               // File containing hydrodynamic data
    std::string type = "hydrodynamic-body";  // Type of the body
    double mass;                             // Mass of the body
    std::vector<double> inertia;             // Inertia of the body
    std::vector<double> position;            // Position of the body
    bool fixed = false;                      // Indicates if the body is fixed
    // For mesh bodies
    std::string geometryFile;  // File path for the geometry
    // For sphere bodies
    double radius;  // Radius of the sphere
};

// Define a struct for PTO configurations
struct PTOConfig {
    std::string type;                              // Type of PTO
    std::vector<int> bodies;                       // Bodies connected to the PTO
    std::vector<std::vector<double>> attachments;  // Attachment points
    double stiffness;                              // Stiffness of the PTO
    double damping;                                // Damping factor

    // Optional fields for specific joint types
    // For linear spring-damper
    double rest_length;            // Rest length of the spring
    double spring_radius  = 0.5;   // Radius of the spring (for visualization)
    int spring_resolution = 1280;  // Resolution of the spring shape
    int spring_turns      = 16;    // Number of turns in the spring shape

    // For rotational PTO
    ChVector<> rotation_axis = ChVector<>(0, 0, 1);  // Default axis
    double initial_angle     = 0;                    // Initial angle (radians)
    std::vector<double> location;                    // Location of the rotational joint
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

    const int colWidth = 24;  // Adjust this value for alignment
    std::cout << "  " << std::left << std::setw(colWidth) << "Model Name:" << config.modelName << "\n";
    std::cout << "  " << std::left << std::setw(colWidth) << "Simulation Duration:" << config.endTime << " seconds\n";
    std::cout << "  " << std::left << std::setw(colWidth) << "Time Step:" << config.dt << " seconds\n\n";

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

    const int colWidth = 24;
    std::cout << "  " << std::left << std::setw(colWidth) << "Type:" << waveConfig.type << "\n";
    std::cout << "  " << std::left << std::setw(colWidth) << "Height:" << waveConfig.height << " meters\n";
    std::cout << "  " << std::left << std::setw(colWidth) << "Period:" << waveConfig.period << " seconds\n";
    if (!waveConfig.spectrumType.empty()) {
        std::cout << "  " << std::left << std::setw(colWidth) << "Spectrum Type:" << waveConfig.spectrumType
                  << std::endl;
    }
    if (!waveConfig.direction.empty()) {
        std::cout << "  " << std::left << std::setw(colWidth) << "Directions:";
        for (const auto& dir : waveConfig.direction) {
            std::cout << dir << " ";
        }
        std::cout << std::endl;
    }
    if (!waveConfig.elevationFile.empty()) {
        std::cout << "  " << std::left << std::setw(colWidth) << "Elevation File:" << waveConfig.elevationFile
                  << std::endl;
    }

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
            // std::cout << "HydroDataFile: " << body.hydroDataFile << std::endl;

        } else if (line.find(".type = '") != std::string::npos) {
            extractAndAssignString(line, body.type, directory, false);
        } else if (line.find(".radius = ") != std::string::npos) {
            body.radius = extractAndAssignDouble(line, ".radius = ");
        } else if (line.find(".geometryFile = '") != std::string::npos) {
            // Extract and assign geometryFile
            extractAndAssignString(line, body.geometryFile, directory);
        } else if (line.find(".mass = ") != std::string::npos) {
            body.mass = extractAndAssignDouble(line, ".mass = ");
        } else if (line.find(".inertia = [") != std::string::npos) {
            extractAndAssignVector(line, body.inertia);
        } else if (line.find(".position = [") != std::string::npos) {
            extractAndAssignVector(line, body.position);
        } else if (line.find(".fixed = '") != std::string::npos) {
            std::string fixedValue = line.substr(line.find("'") + 1, line.rfind("'") - line.find("'") - 1);
            body.fixed             = (fixedValue == "True" || fixedValue == "true");
        }
    }

    // Save the last body
    if (bodyNumber > 0) {
        bodies[bodyNumber] = body;
    }

    // Print summary of bodies parsed
    if (bodies.empty()) {
        std::cout << "No body configurations found in the file." << std::endl;
    } else {
        std::string bodyCountStr = bodies.size() == 1 ? " body" : " bodies";
        std::cout << "  " << bodies.size() << bodyCountStr << " found:" << std::endl;
        for (const auto& [num, config] : bodies) {
            std::string fixedStatus = config.fixed ? "fixed" : "not fixed";
            std::cout << "    - Body " << num << " (" << config.type << ", " << fixedStatus << ")" << std::endl;
        }
    }

    return bodies;
}

// Utility function to parse a vector from a line
std::vector<double> parseVector(const std::string& line) {
    size_t start          = line.find("[") + 1;
    size_t end            = line.find("]");
    std::string vectorStr = line.substr(start, end - start);
    std::istringstream iss(vectorStr);
    std::vector<double> vec;
    double val;
    while (iss >> val) {
        vec.push_back(val);
        if (iss.peek() == ',') iss.ignore();
    }
    return vec;
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
        } else if (line.find(".location = [") != std::string::npos && currentConfig.attachments.empty()) {
            currentConfig.location = parseVector(line);
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
        std::cout << "  " << configs.size() << " PTOs found:" << std::endl;
        for (const auto& config : configs) {
            std::cout << "    - " << config.type << " between body " << config.bodies[0] << " and body "
                      << config.bodies[1] << std::endl;
            if (!config.location.empty()) {
                std::cout << "      Location: [" << config.location[0] << ", " << config.location[1] << ", "
                          << config.location[2] << "]" << std::endl;
            }
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

std::shared_ptr<ChBody> createSphereBody(int bodyNumber,
                                         const BodyConfig& bodyConfig,
                                         ChSystem& system,
                                         const ChColor& color = ChColor(0.244f, 0.225f, 0.072f),  // Default color white
                                         float opacity        = 0.5f) {
    auto sphere = chrono_types::make_shared<ChBodyEasySphere>(bodyConfig.radius,
                                                              1000.0,  // Density
                                                              true,    // Visualization
                                                              false);  // Collisions
    if (!sphere) {
        return nullptr;
    }
    sphere->SetMass(bodyConfig.mass);
    sphere->SetInertiaXX(ChVector<>(bodyConfig.inertia[0], bodyConfig.inertia[1], bodyConfig.inertia[2]));
    sphere->SetPos(ChVector<>(bodyConfig.position[0], bodyConfig.position[1], bodyConfig.position[2]));
    sphere->SetNameString("body" + std::to_string(bodyNumber));

    // Create a visualization material
    auto material = chrono_types::make_shared<ChVisualMaterial>();
    material->SetDiffuseColor(color);
    material->SetOpacity(opacity);
    sphere->GetVisualShape(0)->SetMaterial(0, material);

    system.AddBody(sphere);
    return sphere;
}

std::shared_ptr<ChBody> createSingleBody(int bodyNumber,
                                         const BodyConfig& bodyConfig,
                                         ChSystem& system,
                                         const ChColor& color = ChColor(0.244f, 0.225f, 0.072f),  // Default color
                                         float opacity        = 0.5f) {
    auto body = chrono_types::make_shared<ChBodyEasyMesh>(bodyConfig.geometryFile,
                                                          1000.0,  // Density
                                                          false,   // Automatically evaluate mass
                                                          true,    // Visualization
                                                          false);  // Collisions
    if (!body) {
        std::cerr << "Failed to create body " << bodyNumber << std::endl;
        return nullptr;
    }
    body->SetNameString("body" + std::to_string(bodyNumber));
    setBodyPosition(body, bodyConfig.position);
    body->SetMass(bodyConfig.mass);
    setBodyInertia(body, bodyConfig.inertia);

    // Create a visualization material
    auto material = chrono_types::make_shared<ChVisualMaterial>();
    material->SetDiffuseColor(color);
    material->SetOpacity(opacity);
    body->GetVisualShape(0)->SetMaterial(0, material);

    system.AddBody(body);
    return body;
}

std::pair<std::vector<std::shared_ptr<ChBody>>, std::vector<std::shared_ptr<ChBody>>> createBodies(
    const std::map<int, BodyConfig>& bodyConfigs,
    ChSystem& system) {
    std::vector<std::shared_ptr<ChBody>> allBodies;
    std::vector<std::shared_ptr<ChBody>> hydrodynamicBodies;
    std::shared_ptr<ChBody> ground;

    // First, create hydrodynamic bodies
    for (const auto& [bodyNumber, bodyConfig] : bodyConfigs) {
        if (bodyConfig.type != "nonhydro-sphere") {
            auto body = createSingleBody(bodyNumber, bodyConfig, system);
            if (body) {
                hydrodynamicBodies.push_back(body);
                allBodies.push_back(body);
                std::cout << "    Created hydrodynamic body " << bodyNumber << ".\n";
            }
        }
    }

    // Then, create non-hydrodynamic bodies
    for (const auto& [bodyNumber, bodyConfig] : bodyConfigs) {
        if (bodyConfig.type == "nonhydro-sphere") {
            auto body = createSphereBody(bodyNumber, bodyConfig, system);
            if (body) {
                allBodies.push_back(body);
                std::cout << "    Created non-hydrodynamic body " << bodyNumber << ".\n";
            }
        }
    }

    // Create a ground body if there are any fixed bodies
    bool anyBodyFixed =
        std::any_of(bodyConfigs.begin(), bodyConfigs.end(), [](const auto& pair) { return pair.second.fixed; });

    if (anyBodyFixed) {
        std::cout << "    Creating ground for fixed bodies.\n";
        ground = chrono_types::make_shared<ChBody>();
        system.AddBody(ground);
        ground->SetBodyFixed(true);
        ground->SetCollide(false);
    }

    // Set bodies as fixed if necessary and anchor to the ground
    for (auto& body : allBodies) {
        if (bodyConfigs.at(std::stoi(body->GetNameString().substr(4))).fixed) {
            auto anchor = chrono_types::make_shared<ChLinkMateGeneric>();
            anchor->Initialize(body, ground, false, body->GetFrame_REF_to_abs(), body->GetFrame_REF_to_abs());
            system.Add(anchor);
            anchor->SetConstrainedCoords(true, true, true, true, true, true);
            std::cout << "    Body " << body->GetNameString() << " set as fixed and anchored to ground.\n";
        }
    }

    const int colWidth = 32;
    std::cout << std::setw(colWidth) << std::left << "    Total bodies created:" << allBodies.size() << std::endl;
    std::cout << std::setw(colWidth) << "    Total hydrodynamic bodies:" << hydrodynamicBodies.size() << std::endl;

    return std::make_pair(allBodies, hydrodynamicBodies);
}

void CreateTranslationalPTO(ChSystem& system,
                            const std::vector<std::shared_ptr<ChBody>>& bodies,
                            const PTOConfig& pto_config,
                            std::vector<std::shared_ptr<ChLinkBase>>& ptos) {
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
    prismatic_pto->AddVisualShape(chrono_types::make_shared<ChSpringShape>(
        pto_config.spring_radius, pto_config.spring_resolution, pto_config.spring_turns));
    system.AddLink(prismatic_pto);
    ptos.push_back(prismatic_pto);

    const int colWidth = 24;
    std::cout << "    Created TranslationalPTO between bodies " << pto_config.bodies[0] << " and "
              << pto_config.bodies[1] << std::endl
              << std::setw(colWidth) << std::left << "        Attachments:"
              << "[" << pto_config.attachments[0][0] << ", " << pto_config.attachments[0][1] << ", "
              << pto_config.attachments[0][2] << "] and [" << pto_config.attachments[1][0] << ", "
              << pto_config.attachments[1][1] << ", " << pto_config.attachments[1][2] << "]" << std::endl
              << std::setw(colWidth) << "        Stiffness:" << pto_config.stiffness << std::endl
              << std::setw(colWidth) << "        Damping:" << pto_config.damping << std::endl;
}

void CreateLinSpringDamper(ChSystem& system,
                           const std::vector<std::shared_ptr<ChBody>>& bodies,
                           const PTOConfig& pto_config,
                           std::vector<std::shared_ptr<ChLinkBase>>& ptos) {
    auto body1 = bodies[pto_config.bodies[0] - 1];
    auto body2 = bodies[pto_config.bodies[1] - 1];
    ChVector<> attachment1(pto_config.attachments[0][0], pto_config.attachments[0][1], pto_config.attachments[0][2]);
    ChVector<> attachment2(pto_config.attachments[1][0], pto_config.attachments[1][1], pto_config.attachments[1][2]);

    auto spring_damper = chrono_types::make_shared<ChLinkTSDA>();
    spring_damper->Initialize(body1, body2, true, attachment1, attachment2);
    spring_damper->SetSpringCoefficient(pto_config.stiffness);
    spring_damper->SetDampingCoefficient(pto_config.damping);
    spring_damper->AddVisualShape(chrono_types::make_shared<ChSpringShape>(
        pto_config.spring_radius, pto_config.spring_resolution, pto_config.spring_turns));
    system.AddLink(spring_damper);
    ptos.push_back(spring_damper);

    const int colWidth = 24;
    std::cout << "    Created LinSpringDamper between bodies " << pto_config.bodies[0] << " and "
              << pto_config.bodies[1] << std::endl
              << std::setw(colWidth) << std::left << "        Attachments:"
              << "[" << pto_config.attachments[0][0] << ", " << pto_config.attachments[0][1] << ", "
              << pto_config.attachments[0][2] << "] and [" << pto_config.attachments[1][0] << ", "
              << pto_config.attachments[1][1] << ", " << pto_config.attachments[1][2] << "]" << std::endl
              << std::setw(colWidth) << "        Rest length:" << pto_config.rest_length << std::endl
              << std::setw(colWidth) << "        Stiffness:" << pto_config.stiffness << std::endl
              << std::setw(colWidth) << "        Damping:" << pto_config.damping << std::endl;
}

void CreateRotationalPTO(ChSystem& system,
                         const std::vector<std::shared_ptr<ChBody>>& bodies,
                         const PTOConfig& pto_config,
                         std::vector<std::shared_ptr<ChLinkBase>>& ptos) {
    if (pto_config.bodies.size() != 2) {
        std::cerr << "Error: RotationalPTO requires exactly two bodies." << std::endl;
        return;
    }

    auto body1 = bodies[pto_config.bodies[0] - 1];
    auto body2 = bodies[pto_config.bodies[1] - 1];

    if (pto_config.location.size() != 3) {
        std::cerr << "Error: RotationalPTO 'location' vector must have exactly 3 elements." << std::endl;
        return;
    }

    ChVector<> attachmentPoint(pto_config.location[0], pto_config.location[1], pto_config.location[2]);
    ChQuaternion<> revoluteRot = Q_from_AngX(CH_C_PI / 2.0);

    auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
    revolute->Initialize(body1, body2, ChCoordsys<>(attachmentPoint, revoluteRot));
    system.AddLink(revolute);
    ptos.push_back(revolute);

    const int colWidth = 24;
    std::cout << "    Created RotationalPTO between bodies " << pto_config.bodies[0] << " and " << pto_config.bodies[1]
              << std::endl
              << std::setw(colWidth) << "        Location: "
              << "[" << pto_config.location[0] << ", " << pto_config.location[1] << ", " << pto_config.location[2]
              << "]" << std::endl
              << std::setw(colWidth) << "        Stiffness: " << pto_config.stiffness << std::endl
              << std::setw(colWidth) << "        Damping: " << pto_config.damping << std::endl;
}

void CreateJointOrPTO(ChSystem& system,
                      const std::vector<std::shared_ptr<ChBody>>& bodies,
                      const PTOConfig& pto_config,
                      std::vector<std::shared_ptr<ChLinkBase>>& ptos) {
    if (pto_config.type == "TranslationalPTO") {
        CreateTranslationalPTO(system, bodies, pto_config, ptos);
    } else if (pto_config.type == "LinSpringDamper") {
        CreateLinSpringDamper(system, bodies, pto_config, ptos);
    } else if (pto_config.type == "RotationalPTO") {
        CreateRotationalPTO(system, bodies, pto_config, ptos);
    }
    // Add more else-if clauses for other joint types
}

std::shared_ptr<WaveBase> setupWaveParameters(const WaveConfig& waveConfig,
                                              double timestep,
                                              double simulationDuration,
                                              double rampTime,
                                              size_t num_bodies) {
    std::shared_ptr<WaveBase> hydro_inputs;

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
                       // const std::vector<std::vector<double>>& ptoVelocities,
                       const std::vector<std::vector<double>>& ptoPowers) {
    for (size_t i = 0; i < time_vector.size(); ++i) {
        ptoOutputFile << std::left << std::setw(20) << std::setprecision(2) << std::fixed << time_vector[i];

        // Assuming each inner vector in ptoVelocities and ptoPowers corresponds to a specific PTO
        for (size_t j = 0; j < ptoPowers.size(); ++j) {
            // if (i < ptoVelocities[j].size()) {
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

void printSystemMassMatrix(chrono::ChSystem& system) {
    // Create a sparse matrix to hold the mass matrix
    chrono::ChSparseMatrix massMatrix;

    // Fill the mass matrix
    system.GetMassMatrix(&massMatrix);
    std::cout << "System mass matrix:" << std::endl;
    std::cout << massMatrix << std::endl;

    // Check if the matrix is empty
    if (massMatrix.nonZeros() == 0) {
        std::cout << "Mass matrix is empty." << std::endl;
        return;
    }

    //// Print the mass matrix
    //std::cout << "System Mass Matrix:" << std::endl;
    //for (int k = 0; k < massMatrix.outerSize(); ++k) {
    //    for (Eigen::SparseMatrix<double>::InnerIterator it(massMatrix, k); it; ++it) {
    //        std::cout << "M(" << it.row() << "," << it.col() << ") = " << it.value() << std::endl;
    //    }
    //}
}

// usage: [hydrochrono_wsi.exe] [WECSimInputFile.m] [resultsDirectory] [--gui or --nogui]

int main(int argc, char* argv[]) {
    if (hydroc::SetInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    // Check if the necessary arguments are provided
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <path to .m file> <output directory> [--quiet]" << std::endl;
        return 1;
    }

    // Parse input arguments
    std::string filePath        = argv[1];
    std::string outputDirectory = argv[2];

    // Check for quiet mode argument
    bool quietMode = false;
    std::ofstream logFile;
    std::string logFileName;

    for (int i = 3; i < argc; ++i) {
        if (std::string(argv[i]) == "--quiet") {
            quietMode = true;
            // Create log file name based on input file
            std::filesystem::path inputFilePath(argv[1]);
            logFileName = inputFilePath.parent_path().string();
            if (!logFileName.empty()) {
                logFileName += "/";  // Add a slash if the directory is not empty
            }
            logFileName += inputFilePath.stem().string() + ".log";
            break;
        }
    }

    // Redirect cout to a file if quiet mode is enabled
    if (quietMode) {
        logFile.open(logFileName);
        if (!logFile.is_open()) {
            std::cerr << "Failed to open log file: " << logFileName << std::endl;
            return 1;
        }
        std::cout.rdbuf(logFile.rdbuf());
    }

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

    std::cout << "======== Parsing input file... ========" << std::endl;

    // Parsing simulation configuration
    std::cout << "Simulation Configuration:" << std::endl;
    SimulationConfig simuConfig = parseSimulationConfig(filePath);

    // Parsing wave configuration
    std::cout << "Wave Configuration:" << std::endl;
    WaveConfig waveConfig = parseWaveConfig(filePath);

    // Parsing body information
    std::cout << "\nBody Configuration:" << std::endl;
    std::map<int, BodyConfig> bodyConfigs = parseBodyConfig(filePath);

    // Parsing PTO configuration
    std::cout << "\nPTO Configuration:" << std::endl;
    std::vector<PTOConfig> ptoConfigs = parsePTOConfig(filePath);
    // std::cout << "\n=======================================\n" << std::endl;

    std::cout << "\n\n======== Setting up system...  ========" << std::endl;
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
    auto [bodies, hydrodynamicBodies] = createBodies(bodyConfigs, system);

    std::cout << "\nCreating joint and/or PTO objects..." << std::endl;
    std::vector<std::shared_ptr<ChLinkBase>> ptos;
    for (const auto& ptoConfig : ptoConfigs) {
        if (ptoConfig.type == "RotationalPTO") {
            // Check for either 'location' or 'attachments' for RotationalPTO
            if (!ptoConfig.location.empty()) {
                if (ptoConfig.location.size() != 3) {
                    std::cerr << "Error: RotationalPTO 'location' vector must have exactly 3 elements." << std::endl;
                    continue;  // Skip this iteration if the data is invalid
                }
            } else if (!ptoConfig.attachments.empty()) {
                if (ptoConfig.attachments.size() != 2 || ptoConfig.attachments[0].size() != 3 ||
                    ptoConfig.attachments[1].size() != 3) {
                    std::cerr << "Error: RotationalPTO 'attachments' must each have exactly 3 elements." << std::endl;
                    continue;  // Skip this iteration if the data is invalid
                }
            } else {
                std::cerr << "Error: RotationalPTO requires either 'location' or 'attachments' to be specified."
                          << std::endl;
                continue;  // Skip this iteration if the data is invalid
            }
        } else {
            // For other PTO types, validate the attachments
            if (ptoConfig.attachments.size() != 2 || ptoConfig.attachments[0].size() != 3 ||
                ptoConfig.attachments[1].size() != 3) {
                std::cerr << "Error: PTO attachments must each have exactly 3 elements." << std::endl;
                continue;  // Skip this iteration if the data is invalid
            }
        }

        if (ptoConfig.bodies.size() != 2) {
            std::cerr << "Error: ptoConfig.bodies must contain exactly two body indices." << std::endl;
            continue;  // Skip this iteration if the data is invalid
        }

        // Calling createJointOrPTO() function with the PTOConfig object
        CreateJointOrPTO(system, bodies, ptoConfig, ptos);
    }

    std::cout << "\n\n======== Setting up Waves... ==========\n";
    auto hydro_inputs = setupWaveParameters(waveConfig, simuConfig.dt, simuConfig.endTime, simuConfig.rampTime,
                                            hydrodynamicBodies.size());
    // std::cout << "=======================================\n" << std::endl;

    std::cout << "\n\n======== Initialize HydroChrono... ====\n";

    // std::cerr << "hydroDataFile location: " << bodyConfigs[1].hydroDataFile << std::endl;

    std::cout << "Test hydro..." << std::endl;
    TestHydro hydro_forces(hydrodynamicBodies, bodyConfigs[1].hydroDataFile);
    hydro_forces.AddWaves(hydro_inputs);

    // For profiling
    auto start = std::chrono::high_resolution_clock::now();

    // For output
    std::vector<double> time;
    std::map<int, std::vector<ChVector<>>> body_positions;
    // std::vector<std::vector<double>> ptoVelocities(ptos.size(), std::vector<double>());
    std::vector<std::vector<double>> pto_power(ptos.size(), std::vector<double>());

    // Create an ofstream object
    std::ofstream bodyOutputFile;
    std::ofstream ptoOutputFile;

    // Initialize the output file
    std::cout << "Initialize output files..." << std::endl;
    initializeBodyOutputFile(bodyOutputFile, bodies, outputDirectory);
    initializePTOOutputFile(ptoOutputFile, ptoConfigs.size(), outputDirectory);

    // For visualization

    bool visualizationOn = (simuConfig.explorer == "on");
    // Check for command line arguments to override
    if (argc > 3) {
        if (std::string("--nogui").compare(argv[3]) == 0) {
            visualizationOn = false;
        } else if (std::string("--gui").compare(argv[3]) == 0) {
            visualizationOn = true;
            std::cout << std::endl;
            std::cout << "\n======== Setting up GUI... ============" << std::endl;
            std::cout << std::endl;
        }
    }

    std::cout << "Visualization: " << (visualizationOn ? "Enabled" : "Disabled") << "\n" << std::endl;

    std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);
    hydroc::gui::UI& ui                  = *pui.get();
    ui.Init(&system, simuConfig.modelName.c_str());
    ui.SetCamera(0, -50, -10, 0, 0, -10);
    ui.ShowFrames(true);

    // Main simulation loop
    std::cout << std::endl;
    std::cout << "======== Running simulation... ========" << std::endl;
    std::cout << std::endl;

    // std::vector<std::shared_ptr<ChLinkTSDA>> tsdaPtos;
    // std::vector<std::shared_ptr<ChLinkLockRevolute>> revolutePtos;

    // for (const auto& pto : ptos) {
    //    if (auto tsdaLink = std::dynamic_pointer_cast<ChLinkTSDA>(pto)) {
    //        tsdaPtos.push_back(tsdaLink);
    //    } else if (auto revoluteLink = std::dynamic_pointer_cast<ChLinkLockRevolute>(pto)) {
    //        revolutePtos.push_back(revoluteLink);
    //    }
    //    // Handle other PTO types as needed
    //}

    std::cout << "Number of bodies in system: " << system.Get_bodylist().size() << std::endl;
    
    for (auto& body : system.Get_bodylist()) {
        auto inertia = body->GetInertiaXX();
        std::cout << "Body: " << body->GetNameString() << ", Mass: " << body->GetMass() << ", Inertia: (" << inertia.x()
                  << ", " << inertia.y() << ", " << inertia.z() << ")" << std::endl;
    }

    printSystemMassMatrix(system);

    while (system.GetChTime() <= simulationDuration) {
        if (!ui.IsRunning(timestep)) break;

        if (ui.simulationStarted) {
            system.DoStepDynamics(timestep);
            time.push_back(system.GetChTime());
            // Collect body position data
            for (size_t i = 0; i < bodies.size(); ++i) {
                //body_positions[i].push_back(bodies[i]->GetPos());
            }
            // Collect PTO data for tsdaPtos
            for (size_t i = 0; i < ptos.size(); ++i) {
                // double velocity                = tsdaPtos[i]->GetVelocity();
                // double pto_damping_coefficient = tsdaPtos[i]->GetDampingCoefficient();
                // double power = velocity * velocity * pto_damping_coefficient;  // Assuming pto_damping is accessible
                // ptoVelocities[i].push_back(velocity);
                // ptoPowers[i].push_back(power);
                // std::cout << "==== Running 2... ====" << std::endl;
            }
            // Add similar code for other PTO types if they have different properties or methods
        }
    }

    saveBodyDataToFile(bodyOutputFile, time, body_positions);
    savePTODataToFile(ptoOutputFile, time, pto_power);

    // Restore original cout buffer and close file if it was opened
    if (quietMode) {
        std::cout.rdbuf(std::cout.rdbuf());
        logFile.close();
    }

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