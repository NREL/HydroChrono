#include <hydroc/h5fileinfo.h>
#include <hydroc/helper.h>

#include <cstdlib>
#include <filesystem>  // C++17
#include <iostream>
#include <vector>

using std::filesystem::path;

int main(int argc, char* argv[]) {
    if (hydroc::setInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    path DATADIR(hydroc::getDataDir());

    auto h5fname = (DATADIR / "rm3" / "hydroData" / "rm3.h5").lexically_normal().generic_string();

    HydroData infos = H5FileInfo(h5fname, 2).readH5Data();

    auto rirf_time_vector = infos.GetRIRFTimeVector();

    HydroData infos2 = infos;  // Use move assignement operator

    /* for(auto time: rirf_time_vector) {
         std::cout << time << "\n";
     } */

    std::cout << "End" << std::endl;
    return 0;
}