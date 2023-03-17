#include <hydroc/helper.h>

#include <cstdlib>
#include <filesystem>  // C++17
#include <iostream>
#include <memory>
#include <vector>

using std::filesystem::path;

static path DATADIR{};

int hydroc::setInitialEnvironment(int argc, char* argv[]) noexcept {
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