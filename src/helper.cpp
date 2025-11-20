#include <hydroc/config.h>
#include <hydroc/helper.h>

#include <chrono/core/ChDataPath.h>

#include <cstdlib>
#include <memory>
#include <vector>

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

int hydroc::SetInitialEnvironment(int argc, char* argv[]) noexcept {
    const char* env_p = std::getenv("HYDROCHRONO_DATA_DIR");

    if (env_p) {
        // Highest priority: explicit environment override
        DATADIR = absolute(path(env_p));
        hydroc::cli::LogInfo(std::string("HYDROCHRONO_DATA_DIR set, using '") + getDataDir() + "'");
    } else {
        // If first positional argument looks like a path (does not start with '-'),
        // treat it as the data directory. Otherwise, fall back to the compiled-in
        // default HC_DATA_DIR (build/install data tree).
        if (argc >= 2 && argv[1] && argv[1][0] != '-') {
            DATADIR = absolute(path(argv[1]));
            hydroc::cli::LogInfo(std::string("Using data directory from CLI: '") + getDataDir() + "'");
        } else {
            DATADIR = absolute(path(HC_DATA_DIR));
            hydroc::cli::LogInfo(std::string("Using default data directory HC_DATA_DIR='") + getDataDir() + "'");
        }
    }

    // Set Chrono data directory
    if (std::filesystem::exists(path(std::string(HC_DATA_DIR) + "/chrono"))) {
        chrono::SetChronoDataPath(std::string(HC_DATA_DIR) + "/chrono/");
    } else {
        chrono::SetChronoDataPath(CHRONO_DATA_DIR);
    }

    return 0;
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
