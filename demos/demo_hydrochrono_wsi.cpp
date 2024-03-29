#include <hydroc/gui/guihelper.h>
#include <hydroc/helper.h>
#include <hydroc/hydro_forces.h>

#include <chrono/core/ChRealtimeStep.h>
#include <chrono/physics/ChLinkMate.h>  // Fixed body uses link

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace chrono;
using namespace chrono::geometry;

//////// Parsing data structures and functions

// Struct for simulation configuration
struct SimulationConfig {
    std::string explorer_;    // Explorer type
    std::string model_name_;  // Model name
    double start_time_;       // Start time
    double ramp_time_;        // Ramp time
    double end_time_;         // End time
    double time_step_;        // Time step

    // Camera parameters
    double camera_x_     = 0;
    double camera_y_     = -50;
    double camera_z_     = 5;
    double camera_pitch_ = 0;
    double camera_yaw_   = 0;
    double camera_roll_  = 0;
};


// Struct for wave configuration
struct WaveConfig {
    std::string type_;               // Type of wave
    double height_ = 0.0;            // Height of the wave
    double period_ = 0.0;            // Period of the wave
    std::string spectrum_type_;      // Type of spectrum used
    std::vector<double> direction_;  // Direction of the wave
    std::string elevation_file_;     // File containing elevation data
};

// Struct for body configurations
struct BodyConfig {
    std::string hydro_data_file_;                // File containing hydrodynamic data
    std::string type_    = "hydrodynamic-body";  // Type of the body
    bool is_hydrodynamic = true;                 // Is it a hydrodynamic body or not?
    double mass_;                                // Mass of the body
    std::vector<double> inertia_;                // Inertia of the body
    std::vector<double> position_;               // Position of the body
    bool is_fixed_ = false;                      // Indicates if the body is fixed
    // For mesh bodies
    std::string geometry_file_;  // relpath_;  // Relative file path for the geometry
    // std::string geometry_file_abspath_;  // Absolute file path for the geometry
    // For sphere bodies
    double radius_;  // Radius of the sphere
};

// Struct for PTO (Power Take-Off) configurations
struct PTOConfig {
    std::string type_;                                    // Type of PTO
    std::vector<int> connected_bodies_;                   // Bodies connected to the PTO
    std::vector<std::vector<double>> attachment_points_;  // Attachment points
    double stiffness_;                                    // Stiffness of the PTO
    double damping_;                                      // Damping factor
    std::vector<double> location_;                        // Location of the joint

    // Optional fields for specific joint types
    // For linear spring-damper
    std::optional<double> rest_length_;  // Rest length of the spring
    double spring_radius_  = 0.05;        // Radius of the spring (for visualization)
    int spring_resolution_ = 1280;       // Resolution of the spring shape
    int spring_turns_      = 16;         // Number of turns in the spring shape
    std::optional<double> pretension_;   // Pre-tension in the spring

    // For rotational PTO
    ChVector<> rotation_axis_ = ChVector<>(0, 0, 1);  // Default rotation axis
    double initial_angle_     = 0;                    // Initial angle (radians)
};

// Struct for collating the parsed data structs
struct AllConfigurations {
    SimulationConfig simu_config;
    WaveConfig wave_config;
    std::map<int, BodyConfig> body_configs;
    std::vector<PTOConfig> pto_configs;
};

//// Parser - helper functions

// Helper function to extract a string within single quotes and conditionally prepend a directory
std::string ExtractAndAssignString(const std::string& line,
                                   const std::string& directory = "",
                                   bool append_directory        = true) {
    size_t start          = line.find("'") + 1;
    size_t end            = line.find_last_of("'");
    std::string extracted = line.substr(start, end - start);

    // Prepend the directory only if the extracted path is relative and appending is enabled
    if (append_directory && !directory.empty() && !std::filesystem::path(extracted).is_absolute()) {
        extracted = directory + "/" + extracted;
    }
    return extracted;
}

// Helper function to extract and assign double values
double ExtractAndAssignDouble(const std::string& line, const std::string& pattern) {
    size_t start = line.find(pattern) + pattern.length();
    return std::stod(line.substr(start));
}

// Helper function to extract and return a vector of doubles from a string
std::vector<double> ExtractAndAssignVector(const std::string& line) {
    size_t start        = line.find("[") + 1;
    size_t end          = line.find("]");
    std::string numbers = line.substr(start, end - start);

    // Replacing commas with spaces to handle both as delimiters
    std::replace(numbers.begin(), numbers.end(), ',', ' ');

    std::istringstream iss(numbers);
    std::vector<double> vector;
    double val;
    while (iss >> val) {
        vector.push_back(val);
    }
    return vector;
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

// Helper function to extract a string between single quotes
std::string ExtractStringBetweenQuotes(const std::string& line) {
    size_t start = line.find("'") + 1;
    size_t end   = line.find_last_of("'");
    return line.substr(start, end - start);
}

// Helper function to extract bodies from a string
std::vector<int> ExtractBodies(const std::string& line) {
    std::vector<int> connected_bodies;
    size_t start = line.find("[") + 1;
    size_t end   = line.find("]");

    if (start != std::string::npos && end != std::string::npos) {
        std::string bodies_str = line.substr(start, end - start);
        size_t pos;

        while ((pos = bodies_str.find("body(")) != std::string::npos) {
            size_t start_pos = pos + 5;  // position after "body("
            size_t end_pos   = bodies_str.find(")", start_pos);

            if (end_pos != std::string::npos) {
                std::string body_num_str = bodies_str.substr(start_pos, end_pos - start_pos);
                int body_num             = std::stoi(body_num_str);
                connected_bodies.push_back(body_num);
                bodies_str = bodies_str.substr(end_pos + 1);  // get the remaining string
            } else {
                std::cerr << "Mistake on line in input file defining bodies to be constrained by PTO." << std::endl;
                break;
            }
        }
    }

    return connected_bodies;
}

// Helper function to extract attachment points from a string
std::vector<std::vector<double>> ExtractAttachments(const std::string& line) {
    std::vector<std::vector<double>> attachment_points;
    size_t start                = line.find("[[") + 2;
    size_t end                  = line.find("]]");
    std::string attachments_str = line.substr(start, end - start);

    std::istringstream iss(attachments_str);
    std::string attachment_str;
    while (std::getline(iss, attachment_str, ']')) {
        attachment_str.erase(std::remove(attachment_str.begin(), attachment_str.end(), '['), attachment_str.end());
        attachment_str.erase(std::remove(attachment_str.begin(), attachment_str.end(), ','), attachment_str.end());

        std::istringstream attach_stream(attachment_str);
        std::vector<double> attachment;
        double val;
        while (attach_stream >> val) {
            attachment.push_back(val);
        }
        attachment_points.push_back(attachment);
    }

    return attachment_points;
}

//// Parser - print summary functions

// Helper function to print the summary of parsed body configurations
void PrintParsedBodySummary(const std::map<int, BodyConfig>& body_configs) {
    if (body_configs.empty()) {
        std::cout << "No body configurations found in the file." << std::endl;
        return;
    }

    std::cout << "  " << body_configs.size() << " bodies found:" << std::endl;
    for (const auto& [body_number, config] : body_configs) {
        std::cout << "    - Body " << body_number << " (" << config.type_ << ")";
        if (config.is_fixed_) {
            std::cout << " [Fixed]";
        }
        std::cout << std::endl;

        const int col_width = 32;  // Set the column width

        std::cout << std::setw(col_width) << std::left << "      Mass:" << config.mass_ << std::endl;
        std::cout << std::setw(col_width) << "      Inertia:"
                  << "[" << config.inertia_[0] << ", " << config.inertia_[1] << ", " << config.inertia_[2] << "]"
                  << std::endl;
        std::cout << std::setw(col_width) << "      Position:"
                  << "[" << config.position_[0] << ", " << config.position_[1] << ", " << config.position_[2] << "]"
                  << std::endl;

        if (config.is_hydrodynamic) {
            std::cout << std::setw(col_width) << "      Hydrodynamic Data File:" << config.hydro_data_file_
                      << std::endl;
        }

        if (!config.geometry_file_.empty()) {
            std::cout << std::setw(col_width) << "      Geometry File:" << config.geometry_file_ << std::endl;
        }
    }
}

// Helper function to print the summary of parsed PTO configurations
void PrintParsedPTOSummary(const std::vector<PTOConfig>& configs) {
    if (configs.empty()) {
        std::cout << "No PTO configurations found in the file." << std::endl;
        return;
    }

    std::cout << "  " << configs.size() << " PTOs found:" << std::endl;
    for (const auto& config : configs) {
        std::cout << "    - " << config.type_ << " between body " << config.connected_bodies_[0] << " and body "
                  << config.connected_bodies_[1] << std::endl;
        if (!config.location_.empty()) {
            std::cout << "      Location: [" << config.location_[0] << ", " << config.location_[1] << ", "
                      << config.location_[2] << "]" << std::endl;
        }
    }
}

//// Main parser functions

// Function to parse the Simulation configuration from a given file
SimulationConfig ParseSimulationConfig(const std::string& file_name) {
    // Open the file for reading
    std::ifstream input_file(file_name);
    if (!input_file) {
        // Throw an error if the file can't be opened
        throw std::runtime_error("Failed to open file: " + file_name);
    }

    SimulationConfig config;
    std::string line;
    while (std::getline(input_file, line)) {  // Read each line from the file
        std::istringstream iss(line);
        std::string token;

        // Parse and assign values to the corresponding attributes of the SimulationConfig struct
        if (line.find("simu.modelName = '") != std::string::npos) {
            config.model_name_ = line.substr(line.find("'") + 1, line.rfind("'") - line.find("'") - 1);
        } else if (line.find("simu.explorer = '") != std::string::npos) {
            config.explorer_ = line.substr(line.find("'") + 1, line.rfind("'") - line.find("'") - 1);
        } else if (line.find("simu.startTime =") != std::string::npos) {
            iss >> token >> token >> config.start_time_;
        } else if (line.find("simu.rampTime =") != std::string::npos) {
            iss >> token >> token >> config.ramp_time_;
        } else if (line.find("simu.endTime =") != std::string::npos) {
            iss >> token >> token >> config.end_time_;
        } else if (line.find("simu.dt =") != std::string::npos) {
            iss >> token >> token >> config.time_step_;
        } else if (line.find("simu.cameraX =") != std::string::npos) {
            iss >> token >> token >> config.camera_x_;
        } else if (line.find("simu.cameraY =") != std::string::npos) {
            iss >> token >> token >> config.camera_y_;
        } else if (line.find("simu.cameraZ =") != std::string::npos) {
            iss >> token >> token >> config.camera_z_;
        } else if (line.find("simu.cameraRoll =") != std::string::npos) {
            iss >> token >> token >> config.camera_roll_;
        } else if (line.find("simu.cameraPitch =") != std::string::npos) {
            iss >> token >> token >> config.camera_pitch_;
        } else if (line.find("simu.cameraYaw =") != std::string::npos) {
            iss >> token >> token >> config.camera_yaw_;
        }
    }

    const int kColWidth = 24;  // Adjust this value for alignment
    std::cout << "  " << std::left << std::setw(kColWidth) << "Model Name:" << config.model_name_ << "\n";
    std::cout << "  " << std::left << std::setw(kColWidth) << "Simulation Duration:" << config.end_time_
              << " seconds\n";
    std::cout << "  " << std::left << std::setw(kColWidth) << "Time Step:" << config.time_step_ << " seconds\n\n";

    return config;
}

// Function to parse the Wave configuration from a given file
WaveConfig ParseWaveConfig(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    WaveConfig wave_config;

    // Check if the file is successfully opened
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return wave_config;  // Early return to avoid processing an invalid file
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string word;
        iss >> word;

        // Parse and assign values to the corresponding attributes of the WaveConfig struct
        if (word == "waves") {
            iss >> word;  // Skip '='
            std::string temp;
            iss >> temp;  // Getting the full string, e.g., waveClass('regular');
            // Extracting the word within single quotes
            size_t start      = temp.find("'") + 1;
            size_t end        = temp.find_last_of("'");
            wave_config.type_ = temp.substr(start, end - start);
        } else if (word == "waves.height") {
            iss >> word;  // Skip '='
            iss >> wave_config.height_;
        } else if (word == "waves.period") {
            iss >> word;  // Skip '='
            iss >> wave_config.period_;
        } else if (word == "waves.spectrumType") {
            iss >> word;  // Skip '='
            iss >> wave_config.spectrum_type_;
        } else if (word == "waves.direction") {
            iss >> word;  // Skip '='
            double direction;
            while (iss >> direction) {
                wave_config.direction_.push_back(direction);
            }
        } else if (word == "waves.elevationFile") {
            iss >> word;  // Skip '='
            iss >> wave_config.elevation_file_;
        }
    }

    file.close();

    const int kColWidth = 24;
    std::cout << "  " << std::left << std::setw(kColWidth) << "Type:" << wave_config.type_ << "\n";
    std::cout << "  " << std::left << std::setw(kColWidth) << "Height:" << wave_config.height_ << " meters\n";
    std::cout << "  " << std::left << std::setw(kColWidth) << "Period:" << wave_config.period_ << " seconds\n";
    if (!wave_config.spectrum_type_.empty()) {
        std::cout << "  " << std::left << std::setw(kColWidth) << "Spectrum Type:" << wave_config.spectrum_type_
                  << std::endl;
    }
    if (!wave_config.direction_.empty()) {
        std::cout << "  " << std::left << std::setw(kColWidth) << "Directions:";
        for (const auto& dir : wave_config.direction_) {
            std::cout << dir << " ";
        }
        std::cout << std::endl;
    }
    if (!wave_config.elevation_file_.empty()) {
        std::cout << "  " << std::left << std::setw(kColWidth) << "Elevation File:" << wave_config.elevation_file_
                  << std::endl;
    }

    return wave_config;
}

// Function to parse a Body configuration from a given file
std::map<int, BodyConfig> ParseBodyConfig(const std::string& file_name) {
    std::ifstream input_file(file_name);
    if (!input_file) {
        throw std::runtime_error("Failed to open file.");
    }

    // Extracting directory from the input file path
    std::filesystem::path input_path(file_name);
    std::string directory = input_path.parent_path().string();

    std::map<int, BodyConfig> bodies;
    std::string line;
    BodyConfig body;
    int body_number = 0;  // Start at 0, increment when a new body is found

    while (std::getline(input_file, line)) {
        if (line.find("bodyClass('") != std::string::npos) {
            if (body_number > 0) {  // Save previous body if exists
                bodies[body_number] = body;
            }
            body_number++;        // Increment for new body
            body = BodyConfig();  // Reset body

            // Extract and assign hydroDataFile
            body.hydro_data_file_ = ExtractAndAssignString(line, directory);

        } else if (line.find(".type = '") != std::string::npos) {
            body.type_ = ExtractAndAssignString(line, "", false);
            // Automatically determine if the body is hydrodynamic based on its type
            if (body.type_.find("nonhydro-") == 0) {
                body.is_hydrodynamic = false;
            }
        } else if (line.find(".isHydrodynamic = '") != std::string::npos) {
            std::string hydrodynamic_value = ExtractAndAssignString(line, "", false);
            body.is_hydrodynamic           = (hydrodynamic_value == "True" || hydrodynamic_value == "true");
        } else if (line.find(".radius = ") != std::string::npos) {
            body.radius_ = ExtractAndAssignDouble(line, ".radius = ");
        } else if (line.find(".geometryFile = '") != std::string::npos) {
            body.geometry_file_ = ExtractAndAssignString(line, directory);
        } else if (line.find(".mass = ") != std::string::npos) {
            body.mass_ = ExtractAndAssignDouble(line, ".mass = ");
        } else if (line.find(".inertia = [") != std::string::npos) {
            body.inertia_ = ExtractAndAssignVector(line);
        } else if (line.find(".position = [") != std::string::npos) {
            body.position_ = ExtractAndAssignVector(line);
        } else if (line.find(".fixed = '") != std::string::npos) {
            std::string fixed_value = line.substr(line.find("'") + 1, line.rfind("'") - line.find("'") - 1);
            body.is_fixed_          = (fixed_value == "True" || fixed_value == "true");
        }
    }

    // Save the last body
    if (body_number > 0) {
        bodies[body_number] = body;
    }

    // Print summary of bodies parsed
    if (bodies.empty()) {
        std::cout << "No body configurations found in the file." << std::endl;
    } else {
        PrintParsedBodySummary(bodies);
    }

    return bodies;
}

// Function to parse a PTO configuration from a given file
std::vector<PTOConfig> ParsePTOConfig(const std::string& file_name) {
    std::ifstream input_file(file_name);
    if (!input_file) {
        throw std::runtime_error("Failed to open file.");
    }

    std::vector<PTOConfig> configs;
    PTOConfig current_config;
    std::string line;

    while (std::getline(input_file, line)) {
        if (line.find("ptoClass('") != std::string::npos) {
            // Save the previous config if it is not empty
            if (!current_config.type_.empty()) {
                configs.push_back(current_config);
                current_config = PTOConfig();  // Reset current config
            }
            current_config.type_ = ExtractStringBetweenQuotes(line);
        } else if (line.find(".stiffness = ") != std::string::npos) {
            current_config.stiffness_ = ExtractAndAssignDouble(line, ".stiffness = ");
        } else if (line.find(".damping = ") != std::string::npos) {
            current_config.damping_ = ExtractAndAssignDouble(line, ".damping = ");
        } else if (line.find(".rest_length = ") != std::string::npos) {
            current_config.rest_length_ = ExtractAndAssignDouble(line, ".rest_length = ");
        } else if (line.find(".pretension = ") != std::string::npos) {
            current_config.pretension_ = ExtractAndAssignDouble(line, ".pretension = ");
        } else if (line.find(".bodies = [") != std::string::npos) {
            current_config.connected_bodies_ = ExtractBodies(line);
        } else if (line.find(".attachments = ") != std::string::npos) {
            current_config.attachment_points_ = ExtractAttachments(line);
            // Visualization parameters
        } else if (line.find(".spring_radius = ") != std::string::npos) {
            current_config.spring_radius_ = ExtractAndAssignDouble(line, ".spring_radius = ");
        } else if (line.find(".spring_resolution = ") != std::string::npos) {
            current_config.spring_resolution_ = static_cast<int>(ExtractAndAssignDouble(line, ".spring_resolution = "));
        } else if (line.find(".spring_turns = ") != std::string::npos) {
            current_config.spring_turns_ = static_cast<int>(ExtractAndAssignDouble(line, ".spring_turns = "));
        }
    }
    // Save the last parsed config
    if (!current_config.type_.empty()) {
        configs.push_back(current_config);
    }
    PrintParsedPTOSummary(configs);
    return configs;
}


// Function to parse the input file
AllConfigurations ParseConfigurations(const std::string& file_path) {
    AllConfigurations configs;

    std::cout << "======== Parsing input file... ========" << std::endl;

    std::cout << "Simulation Configuration:" << std::endl;
    configs.simu_config = ParseSimulationConfig(file_path);

    std::cout << "Wave Configuration:" << std::endl;
    configs.wave_config = ParseWaveConfig(file_path);

    std::cout << "\nBody Configuration:" << std::endl;
    configs.body_configs = ParseBodyConfig(file_path);

    std::cout << "\nPTO Configuration:" << std::endl;
    configs.pto_configs = ParsePTOConfig(file_path);

    return configs;
}

// Struct for command line arguments
struct CommandLineArgs {
    std::string file_path;
    std::string output_directory;
    bool quiet_mode;
    std::string log_file_name;
    bool hydro_enabled      = true;
    bool visualization_mode = true;  // true for GUI, false for no GUI
};

//////// Struct and functions to set up the system

// Struct for system and components
struct SystemAndComponents {
    ChSystemNSC system;
    std::vector<std::shared_ptr<ChBody>> all_bodies;
    std::vector<std::shared_ptr<ChBody>> hydrodynamic_bodies;
    std::vector<std::shared_ptr<ChLinkBase>> ptos;
    std::map<int, BodyConfig> body_config_map;  // Map body numbers to their configurations
    std::unique_ptr<TestHydro> hydro_forces;
};

// Helper function to check if parsed PTO configuration is good
bool ValidPTOConfig(const PTOConfig& pto_config) {
    if (pto_config.type_ == "RotationalPTO") {
        if (!pto_config.location_.empty() && pto_config.location_.size() != 3) {
            std::cerr << "Error: RotationalPTO 'location' vector must have exactly 3 elements." << std::endl;
            return false;
        }
        if (!pto_config.attachment_points_.empty() &&
            (pto_config.attachment_points_.size() != 2 || pto_config.attachment_points_[0].size() != 3 ||
             pto_config.attachment_points_[1].size() != 3)) {
            std::cerr << "Error: RotationalPTO 'attachments' must each have exactly 3 elements." << std::endl;
            return false;
        }
        if (pto_config.location_.empty() && pto_config.attachment_points_.empty()) {
            std::cerr << "Error: RotationalPTO requires either 'location' or 'attachments' to be specified."
                      << std::endl;
            return false;
        }
    } else {
        if (pto_config.attachment_points_.size() != 2 || pto_config.attachment_points_[0].size() != 3 ||
            pto_config.attachment_points_[1].size() != 3) {
            std::cerr << "Error: PTO attachments must each have exactly 3 elements." << std::endl;
            return false;
        }
    }

    if (pto_config.connected_bodies_.size() != 2) {
        std::cerr << "Error: ptoConfig.bodies must contain exactly two body indices." << std::endl;
        return false;
    }
    return true;
}

// Helper function to print the system bodies summary
void PrintSystemBodySummary(const SystemAndComponents& sac) {
    if (sac.all_bodies.empty()) {
        std::cout << "No bodies have been created in the system." << std::endl;
        return;
    }

    const int col_width = 20;
    std::cout << "System Body Summary:" << std::endl;

    for (const auto& body : sac.all_bodies) {
        std::string body_name = body->GetNameString();
        int body_number       = -1;

        // Extract body number from name
        if (body_name.rfind("body", 0) == 0 && body_name.length() > 4) {
            body_number = std::stoi(body_name.substr(4));
        }

        // Find the body configuration in the map
        auto it = sac.body_config_map.find(body_number);
        if (it != sac.body_config_map.end()) {
            const auto& body_config = it->second;
            std::cout << std::setw(col_width) << std::left << "  - Body Name:" << body_name << std::endl;
            std::cout << std::setw(col_width) << "    Type:" << body_config.type_ << std::endl;
            std::cout << std::setw(col_width) << "    Mass:" << body->GetMass() << std::endl;
            auto inertia = body->GetInertiaXX();
            std::cout << std::setw(col_width) << "    Inertia:"
                      << "[" << inertia.x() << ", " << inertia.y() << ", " << inertia.z() << "]" << std::endl;
            auto pos = body->GetPos();
            std::cout << std::setw(col_width) << "    Position:"
                      << "[" << pos.x() << ", " << pos.y() << ", " << pos.z() << "]" << std::endl;
        } else {
            std::cerr << "Configuration not found for body: " << body_name << std::endl;
        }
    }
}

//// Functions to set up the Chrono system

// Set body position
void SetBodyPosition(const std::shared_ptr<ChBody>& body, const std::vector<double>& position) {
    if (position.size() == 3) {
        body->SetPos(ChVector<>(position[0], position[1], position[2]));
    } else {
        throw std::runtime_error("Error: Position vector size is not 3.");
    }
}

// Set body inertia
void SetBodyInertia(const std::shared_ptr<ChBody>& body, const std::vector<double>& inertia) {
    if (inertia.size() == 3) {
        body->SetInertiaXX(ChVector<>(inertia[0], inertia[1], inertia[2]));
    } else {
        throw std::runtime_error("Error: Inertia vector size is not 3.");
    }
}

// Create a body for given configuration
std::shared_ptr<ChBody> CreateSingleBody(
    int body_number,
    const BodyConfig& body_config,
    ChSystem& system,
    const ChColor& color = ChColor(0.244f, 0.225f, 0.072f),  // rewrite to make this a body_config attribute
    float opacity        = 0.5f) {
    std::shared_ptr<ChBody> body;

    if (body_config.type_ == "hydrodynamic-body") {
        // For hydrodynamic bodies, use ChBodyEasyMesh
        //std::string absolute_geometry_file = std::filesystem::absolute(body_config.geometry_file_).string();
        body = chrono_types::make_shared<ChBodyEasyMesh>(body_config.geometry_file_,
                                                         1000.0,  // Density
                                                         false,  // Automatically evaluate mass
                                                         true,    // Visualization
                                                         false);  // Collisions
    } else if (body_config.type_ == "nonhydro-sphere") {
        // For sphere bodies, use ChBodyEasySphere
        body = chrono_types::make_shared<ChBodyEasySphere>(body_config.radius_,
                                                           1000.0,  // Density
                                                           true,    // Visualization
                                                           false);  // Collisions
    }

    if (!body) {
        std::cerr << "Failed to create body " << body_number << std::endl;
        return nullptr;
    }

    body->SetNameString("body" + std::to_string(body_number));
    SetBodyPosition(body, body_config.position_);
    body->SetMass(body_config.mass_);
    SetBodyInertia(body, body_config.inertia_);

    auto material = chrono_types::make_shared<ChVisualMaterial>();
    material->SetDiffuseColor(color);
    material->SetOpacity(opacity);
    body->GetVisualShape(0)->SetMaterial(0, material);

    system.AddBody(body);
    return body;
}

// Create all bodies in the system
void CreateBodies(SystemAndComponents& sac, const std::map<int, BodyConfig>& body_configs) {
    std::shared_ptr<ChBody> ground;

    // Create bodies based on configuration
    for (const auto& [body_number, body_config] : body_configs) {
        auto body = CreateSingleBody(body_number, body_config, sac.system);
        if (body) {
            sac.all_bodies.push_back(body);
            if (body_config.is_hydrodynamic) {
                sac.hydrodynamic_bodies.push_back(body);
            }
            sac.body_config_map[body_number] = body_config;  // Make sure this line is present
            // std::cout << "    Created body " << body_number << " of type " << body_config.type_ << ".\n";
        }
    }

    // Create ground if there are any fixed bodies
    bool any_body_fixed =
        std::any_of(body_configs.begin(), body_configs.end(), [](const auto& pair) { return pair.second.is_fixed_; });

    if (any_body_fixed) {
        std::cout << "    Creating ground for fixed bodies.\n";
        ground = chrono_types::make_shared<ChBody>();
        sac.system.AddBody(ground);
        ground->SetPos(ChVector<>(0, 0, 0));
        ground->SetIdentifier(-1);
        ground->SetBodyFixed(true);
        ground->SetCollide(false);

        // Anchor fixed bodies to the ground
        for (auto& body : sac.all_bodies) {
            if (body_configs.at(std::stoi(body->GetNameString().substr(4))).is_fixed_) {
                auto anchor = chrono_types::make_shared<ChLinkMateGeneric>();
                anchor->Initialize(body, ground, false, body->GetFrame_REF_to_abs(), body->GetFrame_REF_to_abs());
                sac.system.Add(anchor);
                anchor->SetConstrainedCoords(true, true, true, true, true, true);
                std::cout << "    Body " << body->GetNameString() << " set as fixed and anchored to ground.\n"
                          << std::endl;
            }
        }
    }

    const int col_width = 32;
    std::cout << std::setw(col_width) << std::left << "    Total bodies created:" << sac.all_bodies.size() << std::endl;
    std::cout << std::setw(col_width) << "    Total hydrodynamic bodies:" << sac.hydrodynamic_bodies.size()
              << std::endl;

    PrintSystemBodySummary(sac);
}

ChVector<> GetChVectorFromConfig(const std::vector<double>& vecConfig) {
    return ChVector<>(vecConfig[0], vecConfig[1], vecConfig[2]);
}

void ConfigurePTOLink(const PTOConfig& pto_config, std::shared_ptr<ChLinkTSDA>& pto_link) {
    pto_link->SetSpringCoefficient(pto_config.stiffness_);
    pto_link->SetDampingCoefficient(pto_config.damping_);

    // Check if rest_length is provided and set it
    if (pto_config.rest_length_.has_value()) {
        pto_link->SetRestLength(pto_config.rest_length_.value());
    }

    // Check if pretension is provided and set it as actuator force
    if (pto_config.pretension_.has_value()) {
        pto_link->SetActuatorForce(pto_config.pretension_.value());
    }

    // Add visual representation of the spring-damper
    pto_link->AddVisualShape(chrono_types::make_shared<ChSpringShape>(
        pto_config.spring_radius_, pto_config.spring_resolution_, pto_config.spring_turns_));
}

std::shared_ptr<ChLinkBase> CreateTranslationalPTO(ChSystem& system,
                                                   const std::vector<std::shared_ptr<ChBody>>& bodies,
                                                   const PTOConfig& pto_config) {
    auto body1             = bodies.at(pto_config.connected_bodies_[0] - 1);
    auto body2             = bodies.at(pto_config.connected_bodies_[1] - 1);
    ChVector<> attachment1 = GetChVectorFromConfig(pto_config.attachment_points_[0]);
    ChVector<> attachment2 = GetChVectorFromConfig(pto_config.attachment_points_[1]);

    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(body1, body2, false, ChCoordsys<>(attachment1), ChCoordsys<>(attachment2));
    system.AddLink(prismatic);

    auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
    prismatic_pto->Initialize(body1, body2, false, attachment1, attachment2);
    ConfigurePTOLink(pto_config, prismatic_pto);
    system.AddLink(prismatic_pto);
    return prismatic_pto;
}

std::shared_ptr<ChLinkBase> CreateLinSpringDamper(ChSystem& system,
                                                  const std::vector<std::shared_ptr<ChBody>>& bodies,
                                                  const PTOConfig& pto_config) {
    auto body1             = bodies.at(pto_config.connected_bodies_[0] - 1);
    auto body2             = bodies.at(pto_config.connected_bodies_[1] - 1);
    ChVector<> attachment1 = GetChVectorFromConfig(pto_config.attachment_points_[0]);
    ChVector<> attachment2 = GetChVectorFromConfig(pto_config.attachment_points_[1]);

    auto spring_damper = chrono_types::make_shared<ChLinkTSDA>();
    spring_damper->Initialize(body1, body2, true, attachment1, attachment2);
    ConfigurePTOLink(pto_config, spring_damper);
    system.AddLink(spring_damper);
    return spring_damper;
}

std::shared_ptr<ChLinkBase> CreateRotationalPTO(ChSystem& system,
                                                const std::vector<std::shared_ptr<ChBody>>& bodies,
                                                const PTOConfig& pto_config) {
    auto body1                 = bodies.at(pto_config.connected_bodies_[0] - 1);
    auto body2                 = bodies.at(pto_config.connected_bodies_[1] - 1);
    ChVector<> attachmentPoint = GetChVectorFromConfig(pto_config.location_);
    ChQuaternion<> revoluteRot = Q_from_AngX(CH_C_PI / 2.0);

    auto revolute_pto = chrono_types::make_shared<ChLinkLockRevolute>();
    revolute_pto->Initialize(body2, body1, ChCoordsys<>(attachmentPoint, revoluteRot));
    system.AddLink(revolute_pto);
    return revolute_pto;
}

void CreateJointOrPTO(SystemAndComponents& sac, const PTOConfig& pto_config) {
    // Check if there are sufficient details to create the PTO
    bool valid_config = true;
    if (pto_config.type_ == "RotationalPTO") {
        valid_config = !pto_config.location_.empty() && pto_config.connected_bodies_.size() == 2;
    } else {
        valid_config = pto_config.connected_bodies_.size() >= 2 && pto_config.attachment_points_.size() >= 2;
    }

    if (!valid_config) {
        std::cerr << "Insufficient bodies or attachments for " << pto_config.type_ << " configuration." << std::endl;
        return;
    }

    // Common logging before creating the PTO
    constexpr int col_width = 24;
    std::cout << "    Creating " << pto_config.type_ << " between bodies " << pto_config.connected_bodies_[0] << " and "
              << pto_config.connected_bodies_[1] << std::endl;
    if (pto_config.type_ != "RotationalPTO") {
        std::cout << std::setw(col_width) << std::left << "        Attachments:"
                  << "[" << pto_config.attachment_points_[0][0] << ", " << pto_config.attachment_points_[0][1] << ", "
                  << pto_config.attachment_points_[0][2] << "] and [" << pto_config.attachment_points_[1][0] << ", "
                  << pto_config.attachment_points_[1][1] << ", " << pto_config.attachment_points_[1][2] << "]"
                  << std::endl;
    }

    std::cout << std::setw(col_width) << "        Stiffness:" << pto_config.stiffness_ << std::endl
              << std::setw(col_width) << "        Damping:" << pto_config.damping_ << std::endl;

    // Creating the specific PTO
    std::shared_ptr<ChLinkBase> pto_link;
    if (pto_config.type_ == "TranslationalPTO") {
        pto_link = CreateTranslationalPTO(sac.system, sac.all_bodies, pto_config);
    } else if (pto_config.type_ == "LinSpringDamper") {
        pto_link = CreateLinSpringDamper(sac.system, sac.all_bodies, pto_config);
    } else if (pto_config.type_ == "RotationalPTO") {
        pto_link = CreateRotationalPTO(sac.system, sac.all_bodies, pto_config);
    } else {
        std::cerr << "Unknown PTO type: " << pto_config.type_ << std::endl;
        return;
    }

    if (pto_link) {
        sac.ptos.push_back(pto_link);
    }
}

void PrintWaveParams(const IrregularWaveParams& params) {
    std::cout << std::left << std::setw(20) << "Simulation Timestep:" << params.simulation_dt_ << "\n"
              << std::left << std::setw(20) << "Simulation Duration:" << params.simulation_duration_ << "\n"
              << std::left << std::setw(20) << "Ramp Time:" << params.ramp_duration_ << "\n"
              << std::left << std::setw(20) << "Wave Height:" << params.wave_height_ << "\n"
              << std::left << std::setw(20) << "Wave Period:" << params.wave_period_ << "\n";
}

std::shared_ptr<WaveBase> SetupWaveParameters(const WaveConfig& wave_config,
                                              double timestep,
                                              double simulation_duration,
                                              double ramp_time,
                                              size_t num_bodies) {
    std::shared_ptr<WaveBase> wave;

    std::cout << std::left << std::setw(20) << "Wave Type: " << wave_config.type_ << "\n";

    if (wave_config.type_ == "irregular") {
        IrregularWaveParams wave_params;
        wave_params.num_bodies_          = num_bodies;
        wave_params.simulation_dt_       = timestep;
        wave_params.simulation_duration_ = simulation_duration;
        wave_params.ramp_duration_       = ramp_time;
        wave_params.wave_height_         = wave_config.height_;
        wave_params.wave_period_         = wave_config.period_;

        PrintWaveParams(wave_params);

        try {
            wave = std::make_shared<IrregularWaves>(wave_params);
        } catch (const std::exception& e) {
            std::cerr << "Exception in IrregularWaves creation: " << e.what() << '\n';
            return nullptr;
        }
    } else if (wave_config.type_ == "regular") {
        auto regular_wave                     = std::make_shared<RegularWave>(num_bodies);
        regular_wave->regular_wave_amplitude_ = (wave_config.height_ / 2);
        regular_wave->regular_wave_omega_     = (2.0 * CH_C_PI / wave_config.period_);

        std::cout << std::left << std::setw(20) << "Wave Amplitude:" << regular_wave->regular_wave_amplitude_ << "\n";
        std::cout << std::left << std::setw(20) << "Wave Omega:" << regular_wave->regular_wave_omega_ << "\n";

        wave = regular_wave;
    } else if (wave_config.type_ == "still" || wave_config.type_ == "noWaveCIC") {
        wave = std::make_shared<NoWave>(num_bodies);
    } else {
        std::cerr << "Unknown wave type: " << wave_config.type_ << std::endl;
    }

    return wave;
}

bool OpenOutputFile(std::ofstream& output_file, const std::string& file_path) {
    std::filesystem::path output_path(file_path);

    std::filesystem::create_directories(output_path.parent_path());

    output_file.open(output_path);
    return output_file.is_open();
}

void InitializeBodyOutputFile(std::ofstream& body_output_file,
                              const std::vector<std::shared_ptr<ChBody>>& bodies,
                              const std::string& output_directory) {
    constexpr size_t kColumnWidth   = 20;
    constexpr size_t kBodyDataWidth = 16;

    std::filesystem::path output_file = std::filesystem::path(output_directory) / "body_output.txt";

    if (OpenOutputFile(body_output_file, output_file.string())) {
        body_output_file << std::left << std::setw(kColumnWidth) << "Time (s)";
        for (size_t i = 0; i < bodies.size(); ++i) {
            std::string body_index = "Body" + std::to_string(i);
            body_output_file << std::setw(kBodyDataWidth) << (body_index + "_x (m)") << std::setw(kBodyDataWidth)
                             << (body_index + "_y (m)") << std::setw(kBodyDataWidth) << (body_index + "_z (m)");
        }
        body_output_file << std::endl;
    } else {
        std::cerr << "Failed to open the body output file at: " << output_file << std::endl;
    }
}

void InitializePTOOutputFile(std::ofstream& pto_output_file, size_t pto_count, const std::string& output_directory) {
    constexpr size_t kColumnWidth = 20;

    std::filesystem::path output_file = std::filesystem::path(output_directory) / "pto_output.txt";

    if (OpenOutputFile(pto_output_file, output_file.string())) {
        pto_output_file << std::left << std::setw(kColumnWidth) << "Time (s)";
        for (size_t i = 0; i < pto_count; ++i) {
            pto_output_file << std::setw(kColumnWidth) << ("PTO" + std::to_string(i) + "_power (W)");
        }
        pto_output_file << std::endl;
    } else {
        std::cerr << "Failed to open the PTO output file at: " << output_file << std::endl;
    }
}

void CollectPTOData(const std::vector<std::shared_ptr<ChLinkTSDA>>& ptos,
                    std::vector<double>& pto_velocities,
                    std::vector<double>& pto_powers) {
    pto_velocities.clear();
    pto_powers.clear();

    for (const auto& pto : ptos) {
        double pto_velocity   = pto->GetVelocity();
        double pto_damping    = pto->GetDampingCoefficient();
        double power_absorbed = pto_velocity * pto_velocity * pto_damping;

        pto_velocities.push_back(pto_velocity);
        pto_powers.push_back(power_absorbed);
    }
}

CommandLineArgs InitializeAndParseArgs(int argc, char* argv[]) {
    if (argc < 3) {
        throw std::runtime_error("Usage: " + std::string(argv[0]) +
                                 " <path to .m file> <output directory> [--quiet] [--nohydro] [--gui/--nogui]");
    }

    CommandLineArgs args;
    args.file_path        = argv[1];
    args.output_directory = argv[2];
    args.quiet_mode       = false;

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--quiet") {
            args.quiet_mode = true;
            std::filesystem::path input_file_path(argv[1]);
            args.log_file_name = input_file_path.parent_path().string();
            if (!args.log_file_name.empty()) {
                args.log_file_name += "/";
            }
            args.log_file_name += input_file_path.stem().string() + ".log";
        } else if (arg == "--nohydro") {
            args.hydro_enabled = false;
        } else if (arg == "--nogui") {
            args.visualization_mode = false;
        } else if (arg == "--gui") {
            args.visualization_mode = true;
        }
    }

    return args;
}

void SetupLogging(const CommandLineArgs& args) {
    if (args.quiet_mode) {
        std::ofstream log_file(args.log_file_name);
        if (!log_file.is_open()) {
            throw std::runtime_error("Failed to open log file: " + args.log_file_name);
        }
        std::cout.rdbuf(log_file.rdbuf());
    }
    // If not in quiet mode, standard output remains unchanged
}

void PrintBanner() {
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
}

void LogChronoVersion() {
    // Assuming CHRONO_VERSION is a predefined macro or a global constant
    // containing the version of the Chrono library.
    GetLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";
}

SystemAndComponents SetupSystem(const AllConfigurations& configurations, bool hydro_enabled) {
    SystemAndComponents sac;

    std::cout << "\n\n======== Setting up system...  ========" << std::endl;

    sac.system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));
    sac.system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    sac.system.SetSolverType(ChSolver::Type::GMRES);
    sac.system.SetStep(configurations.simu_config.time_step_);

    std::cout << "Creating body objects..." << std::endl;
    CreateBodies(sac, configurations.body_configs);

    if (configurations.pto_configs.empty()) {
        std::cout << "No PTO configurations found. Continuing without PTOs." << std::endl;
    } else {
        std::cout << "\nCreating joint and/or PTO objects..." << std::endl;
        for (const auto& pto_config : configurations.pto_configs) {
            if (ValidPTOConfig(pto_config)) {
                CreateJointOrPTO(sac, pto_config);
            }
        }
    }

    // Find the first hydrodynamic body to get .h5 file
    std::string hydro_data_file;
    for (const auto& [key, body_config] : configurations.body_configs) {
        if (body_config.is_hydrodynamic) {
            hydro_data_file = body_config.hydro_data_file_;
            break;  // Found the first hydrodynamic body, exit the loop
        }
    }

    // Setting up waves and initializing hydrodynamic forces
    if (hydro_enabled && !hydro_data_file.empty()) {
        std::cout << "\n\n======== Setting up Waves... ==========\n";
        auto hydro_inputs = SetupWaveParameters(configurations.wave_config, configurations.simu_config.time_step_,
                                                configurations.simu_config.end_time_,
                                                configurations.simu_config.ramp_time_, sac.hydrodynamic_bodies.size());

        std::cout << "\n\n======== Initialize HydroChrono... ====\n";
        sac.hydro_forces = std::make_unique<TestHydro>(sac.hydrodynamic_bodies, hydro_data_file);
        sac.hydro_forces->AddWaves(hydro_inputs);
    } else if (!hydro_enabled) {
        std::cout << "Hydrodynamic forces disabled." << std::endl;
    } else {
        std::cerr << "Error: No hydrodynamic body found in the configuration." << std::endl;
    }

    return sac;
}

std::shared_ptr<hydroc::gui::UI> SetupVisualization(const CommandLineArgs& args,
                                                    const AllConfigurations& configurations,
                                                    SystemAndComponents& sac) {
    // Set visualization based on simu_config's explorer_ flag
    bool visualization_on = (configurations.simu_config.explorer_ == "on");

    // Override with command line argument if provided
    if (!args.visualization_mode) {
        visualization_on = false;
    }

    std::cout << "Visualization: " << (visualization_on ? "Enabled" : "Disabled") << "\n" << std::endl;

    if (!visualization_on) {
        return nullptr;  // If visualization is disabled, return null
    }

    // Initialize UI if visualization is enabled
    auto ui = hydroc::gui::CreateUI(visualization_on);
    ui->Init(&sac.system, configurations.simu_config.model_name_.c_str());
    ui->SetCamera(configurations.simu_config.camera_x_, configurations.simu_config.camera_y_,
                  configurations.simu_config.camera_z_, configurations.simu_config.camera_pitch_,
                  configurations.simu_config.camera_yaw_, configurations.simu_config.camera_roll_);
    ui->ShowFrames(true);

    return ui;  // Return the created UI object
}

class SimulationOutput {
  public:
    SimulationOutput(const std::vector<std::shared_ptr<ChBody>>& bodies,
                     size_t pto_count,
                     const std::string& output_directory)
        : bodies_(bodies),
          pto_powers_(pto_count, std::vector<double>()),
          body_output_file_path_(std::filesystem::absolute(output_directory + "/body_output.txt").string()),
          pto_output_file_path_(std::filesystem::absolute(output_directory + "/pto_output.txt").string()) {}

    void AppendSimulationData(double time,
                              const std::vector<std::shared_ptr<ChBody>>& bodies,
                              const std::vector<std::shared_ptr<ChLinkBase>>& ptos);

    double CalculatePTOPower(const std::shared_ptr<ChLinkBase>& pto);

    void SaveDataToFile() {
        SaveBodyDataToFile();
        SavePTODataToFile();
        std::cout << "Simulation data saved successfully." << std::endl;
        std::cout << "Body data saved to: " << body_output_file_path_ << std::endl;
        std::cout << "PTO data saved to: " << pto_output_file_path_ << std::endl;
    }

  private:
    std::ofstream body_output_file_;
    std::ofstream pto_output_file_;
    std::vector<double> time_vector_;
    std::map<std::string, std::vector<ChVector<>>> body_positions_;
    std::vector<std::vector<double>> pto_powers_;
    std::vector<std::shared_ptr<ChBody>> bodies_;
    std::string body_output_file_path_;
    std::string pto_output_file_path_;

    void SaveBodyDataToFile() {
        std::filesystem::path dir = std::filesystem::path(body_output_file_path_).parent_path();
        if (!std::filesystem::exists(dir)) {
            std::filesystem::create_directories(dir);
        }
        constexpr size_t kTimeColumnWidth = 20;
        constexpr size_t kDataColumnWidth = 16;
        constexpr int kPrecision          = 4;

        body_output_file_.open(body_output_file_path_);
        if (!body_output_file_.is_open()) {
            std::cerr << "Failed to open body output file at " << body_output_file_path_ << std::endl;
            return;
        }

        // Write headers
        body_output_file_ << std::left << std::setw(kTimeColumnWidth) << "Time (s)";
        for (size_t i = 0; i < bodies_.size(); ++i) {
            std::string body_index = "Body" + std::to_string(i);
            body_output_file_ << std::setw(kDataColumnWidth) << (body_index + "_x (m)") << std::setw(kDataColumnWidth)
                              << (body_index + "_y (m)") << std::setw(kDataColumnWidth) << (body_index + "_z (m)");
        }
        body_output_file_ << std::endl;

        // Write data
        for (size_t i = 0; i < time_vector_.size(); ++i) {
            body_output_file_ << std::left << std::setw(kTimeColumnWidth) << std::setprecision(2) << std::fixed
                              << time_vector_[i];
            for (const auto& [body_name, positions] : body_positions_) {
                const ChVector<>& position = positions[i];
                body_output_file_ << std::setw(kDataColumnWidth) << std::setprecision(kPrecision) << std::fixed
                                  << position.x() << std::setw(kDataColumnWidth) << position.y()
                                  << std::setw(kDataColumnWidth) << position.z();
            }
            body_output_file_ << std::endl;
        }

        body_output_file_.close();
    }

    void SavePTODataToFile() {
        std::filesystem::path dir = std::filesystem::path(pto_output_file_path_).parent_path();
        if (!std::filesystem::exists(dir)) {
            std::filesystem::create_directories(dir);
        }
        constexpr size_t kTimeColumnWidth  = 20;
        constexpr size_t kPowerColumnWidth = 20;
        constexpr int kPrecision           = 4;

        pto_output_file_.open(pto_output_file_path_);
        if (!pto_output_file_.is_open()) {
            std::cerr << "Failed to open PTO output file at " << pto_output_file_path_ << std::endl;
            return;
        }

        // Write headers
        pto_output_file_ << std::left << std::setw(kTimeColumnWidth) << "Time (s)";
        for (size_t i = 0; i < pto_powers_.size(); ++i) {
            pto_output_file_ << std::setw(kPowerColumnWidth) << ("PTO" + std::to_string(i) + "_power (W)");
        }
        pto_output_file_ << std::endl;

        // Write data
        for (size_t i = 0; i < time_vector_.size(); ++i) {
            pto_output_file_ << std::left << std::setw(kTimeColumnWidth) << std::setprecision(2) << std::fixed
                             << time_vector_[i];
            for (size_t j = 0; j < pto_powers_.size(); ++j) {
                if (i < pto_powers_[j].size()) {
                    pto_output_file_ << std::setw(kPowerColumnWidth) << std::setprecision(kPrecision) << std::fixed
                                     << pto_powers_[j][i];
                }
            }
            pto_output_file_ << std::endl;
        }

        pto_output_file_.close();
    }
};

double SimulationOutput::CalculatePTOPower(const std::shared_ptr<ChLinkBase>& pto) {
    // Example calculation for a specific type of PTO
    if (auto tsdaLink = std::dynamic_pointer_cast<ChLinkTSDA>(pto)) {
        double velocity            = tsdaLink->GetVelocity();
        double damping_coefficient = tsdaLink->GetDampingCoefficient();
        return velocity * velocity * damping_coefficient;
    }
    // Handle other PTO types as needed
    return 0.0;  // Default return if PTO type is not handled
}

void SimulationOutput::AppendSimulationData(double time,
                                            const std::vector<std::shared_ptr<ChBody>>& bodies,
                                            const std::vector<std::shared_ptr<ChLinkBase>>& ptos) {
    // Append data for bodies
    time_vector_.push_back(time);
    for (auto& body : bodies) {
        body_positions_[body->GetNameString()].push_back(body->GetPos());
    }

    // Append data for PTOs
    for (size_t i = 0; i < ptos.size(); ++i) {
        double power = CalculatePTOPower(ptos[i]);
        pto_powers_[i].push_back(power);
    }
}

void RunSimulation(const AllConfigurations& configurations,
                   SystemAndComponents& sac,
                   hydroc::gui::UI& ui,
                   SimulationOutput* output) {
    using Clock = std::chrono::high_resolution_clock;
    std::cout << "\n======== Running simulation with GUI... ========\n";

    auto start_time = Clock::now();

    // Run the simulation
    try {
        while (sac.system.GetChTime() <= configurations.simu_config.end_time_) {
            if (!ui.IsRunning(configurations.simu_config.time_step_)) {
                std::cout << "Simulation stopped by user.\n";
                break;
            }

            if (ui.simulationStarted) {
                sac.system.DoStepDynamics(configurations.simu_config.time_step_);
                output->AppendSimulationData(sac.system.GetChTime(), sac.all_bodies, sac.ptos);
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception during simulation with GUI: " << e.what() << '\n';
    }

    auto end_time = Clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "Simulation with GUI completed in " << std::fixed << std::setprecision(3) << duration.count()
              << " seconds.\n";
}

void RunSimulation(const AllConfigurations& configurations, SystemAndComponents& sac, SimulationOutput* output) {
    using Clock = std::chrono::high_resolution_clock;
    std::cout << "\n======== Running headless simulation... ========\n";

    auto start_time = Clock::now();

    // Run the simulation
    try {
        while (sac.system.GetChTime() <= configurations.simu_config.end_time_) {
            sac.system.DoStepDynamics(configurations.simu_config.time_step_);
            output->AppendSimulationData(sac.system.GetChTime(), sac.all_bodies, sac.ptos);
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception during headless simulation: " << e.what() << '\n';
    }

    auto end_time = Clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "Headless simulation completed in " << std::fixed << std::setprecision(3) << duration.count()
              << " seconds.\n";
}
// usage: [hydrochrono.exe] [InputFile.m] [resultsDir] [--gui/--nogui] [--quiet] [--nohydro]
int main(int argc, char* argv[]) {
    hydroc::SetInitialEnvironment(argc, argv);

    try {
        auto args = InitializeAndParseArgs(argc, argv);
        SetupLogging(args);

        PrintBanner();
        LogChronoVersion();

        auto configurations = ParseConfigurations(args.file_path);
        auto sac            = SetupSystem(configurations, args.hydro_enabled);

        SimulationOutput output(sac.all_bodies, sac.ptos.size(), args.output_directory);

        auto ui = SetupVisualization(args, configurations, sac);
        if (ui) {
            RunSimulation(configurations, sac, *ui, &output);
        } else {
            RunSimulation(configurations, sac, &output);
        }

        std::cout << "\n======== Saving results... ============" << std::endl;
        output.SaveDataToFile();

        std::cout << "\nSimulation completed successfully." << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}