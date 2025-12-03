/*********************************************************************
 * @file  helper.cpp
 *
 * @brief Implementation of helper utilities.
 *********************************************************************/

#include <hydroc/config.h>
#include <hydroc/helper.h>

#include <chrono/core/ChDataPath.h>
#include <chrono_thirdparty/cxxopts/ChCLI.h>

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

bool hydroc::GetCLIArguments(int argc,
                             char** argv,
                             const std::string& description,
                             bool& output,
                             bool& profile,
                             bool& gui,
                             std::string& data_dir) {
    chrono::ChCLI cli(argv[0], description);

    cli.AddOption<std::string>("", "data_dir", "HydrChrono data directory", "");
    if (output)
        cli.AddOption<bool>("", "no_output", "Disable generation of simulation output");
    else
        cli.AddOption<bool>("", "output", "Enable generation of simulation output");

    if (profile)
        cli.AddOption<bool>("", "no_profile", "Disable profiling of simulation time");
    else
        cli.AddOption<bool>("", "profile", "Enable profiling of simulation time");

    if (gui)
        cli.AddOption<bool>("", "no_gui", "Disable GUI");
    else
        cli.AddOption<bool>("", "gui", "Enable GUI");

    if (!cli.Parse(argc, argv)) {
        cli.Help();
        return false;
    }

    data_dir = cli.Get("data_dir").as<std::string>();

    if (output)
        output = !cli.GetAsType<bool>("no_output");
    else
        output = cli.GetAsType<bool>("output");

    if (profile)
        profile = !cli.GetAsType<bool>("no_profile");
    else
        profile = cli.GetAsType<bool>("profile");

    if (gui)
        gui = !cli.GetAsType<bool>("no_gui");
    else
        gui = cli.GetAsType<bool>("gui");

    return true;
}

size_t get_lower_index(double value, const std::vector<double>& ticks) {
    auto it = std::upper_bound(ticks.begin(), ticks.end(), value);
    // get nearest-below index
    size_t idx = it - ticks.begin() - 1;
    // remove one if equal to value
    if (ticks[idx] == value) {
        idx -= 1;
    }
    if (idx <= 0 || idx >= ticks.size() - 1) {
        throw std::runtime_error("Could not find index for value " + std::to_string(value) + " in array with bounds (" +
                                 std::to_string(ticks.front()) + ", " + std::to_string(ticks.back()) + ").");
    }
    // return index
    return idx;
}

using std::filesystem::path;

static path DATADIR{};

bool hydroc::SetInitialEnvironment(const std::string& data_dir) noexcept {
    const char* env_p = std::getenv("HYDROCHRONO_DATA_DIR");

    if (env_p) {
        // Highest priority: explicit environment override
        DATADIR = absolute(path(env_p));
        hydroc::cli::LogInfo(std::string("Using data directory from HYDROCHRONO_DATA_DIR: '") + getDataDir() + "'");
    } else if (!data_dir.empty()) {
        DATADIR = absolute(path(data_dir));
        hydroc::cli::LogInfo(std::string("Using provided data directory: '") + getDataDir() + "'");
    } else {
        DATADIR = absolute(path(HC_DATA_DIR));
        hydroc::cli::LogInfo(std::string("Using default data directory HC_DATA_DIR: '") + getDataDir() + "'");
    }

    // Set Chrono data directory
    if (std::filesystem::exists(path(std::string(HC_DATA_DIR) + "/chrono"))) {
        chrono::SetChronoDataPath(std::string(HC_DATA_DIR) + "/chrono/");
    } else {
        chrono::SetChronoDataPath(CHRONO_DATA_DIR);
    }

    return true;
}

std::string hydroc::getDataDir() noexcept {
    return DATADIR.lexically_normal().generic_string();
}

void hydroc::ensure_directory_exists(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        std::cout << "Path " << std::filesystem::absolute(path) << " does not exist, creating it now..." << std::endl;
        std::filesystem::create_directory(path);
    }
}

std::string hydroc::getDemoOutDir() {
    std::string dir_name = "results/demos";
    std::filesystem::create_directories(dir_name);
    return dir_name;
}

std::string hydroc::getTestOutDir() {
    std::string dir_name = "results/tests";
    std::filesystem::create_directories(dir_name);
    return dir_name;
}
