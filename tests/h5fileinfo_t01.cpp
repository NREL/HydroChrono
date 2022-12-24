#include <hydroc/h5fileinfo.h>

#include <filesystem>  // C++17
#include <cstdlib>
#include <iostream>
#include <vector>

using std::filesystem::path;
using std::filesystem::absolute;

static path DATADIR{};

int main(int argc, char* argv[]) {

    const char* env_p = std::getenv("HYDRO_CHRONO_DATA_DIR");

    if (env_p == nullptr) {
        if (argc < 2) {
            std::cerr << "Usage: .exe [<datadir>] or set HYDRO_CHRONO_DATA_DIR environement variable" << std::endl;
            return 1;
        } else {
            DATADIR = absolute(path(argv[1]));
        }
    } else {
        DATADIR = absolute(path(env_p));
    }

    //auto h5fname = (DATADIR / "sphere" / "hydroData" /"sphere.h5").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "rm3" / "hydroData" /"rm3.h5").lexically_normal().generic_string();


    //H5FileInfo(h5fname, "body1");

    std::vector<H5FileInfo> infos;
    infos.emplace_back(h5fname, "body1");
    infos.emplace_back(h5fname, "body2"); // reopen the file

    auto rirf_time_vector = infos[0].GetRIRFTimeVector();

    for(auto time: rirf_time_vector) {
        std::cout << time << "\n";
    }

    std::cout << "End" << std::endl;
    return 0;
}