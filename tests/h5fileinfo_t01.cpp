#include <hydroc/h5fileinfo.h>
#include <hydroc/helper.h>

#include <filesystem>  // C++17
#include <cstdlib>
#include <iostream>
#include <vector>

using std::filesystem::path;

int main(int argc, char* argv[]) {

    if (hydroc::setInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    path DATADIR(hydroc::getDataDir());


    auto h5fname = (DATADIR / "rm3" / "hydroData" /"rm3.h5").lexically_normal().generic_string();


    //H5FileInfo(h5fname, "body1");

    std::vector<H5FileInfo> infos;
    infos.emplace_back(h5fname, "body1"); // Use move constructor
    infos.emplace_back(h5fname, "body2"); // reopen the file

    auto rirf_time_vector = infos[0].GetRIRFTimeVector();

    std::vector<H5FileInfo> infos2 = infos; // Use move assignement operator

   /* for(auto time: rirf_time_vector) {
        std::cout << time << "\n";
    } */

    std::cout << "End" << std::endl;
    return 0;
}