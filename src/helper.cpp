#include <hydroc/helper.h>

#include <cstdlib>
#include <filesystem>  // C++17
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

    if (env_p == nullptr) {
        if (argc < 2) {
            std::cerr << "Warning::Usage: .exe [<datadir>] or set HYDROCHRONO_DATA_DIR environement variable"
                      << std::endl;

            DATADIR = absolute(path("..") / ".." / "demos");
            std::cerr << "Set default demos path to'" << getDataDir() << "'" << std::endl;
            return 0;
        } else {
            DATADIR = absolute(path(argv[1]));
        }
    } else {
        DATADIR = absolute(path(env_p));
    }
    return 0;
}

std::string hydroc::getDataDir() noexcept {
    return DATADIR.lexically_normal().generic_string();
}