#include <hydroc/h5fileinfo.h>
#include <hydroc/chloadaddedmass.h>
#include <hydroc/helper.h>

#include <chrono/core/ChTypes.h>
#include <chrono/physics/ChBodyEasy.h>

#include <filesystem>  // C++17
#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>

using std::filesystem::path;
using std::filesystem::absolute;

int main(int argc, char* argv[]) {
    if (hydroc::setInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    path DATADIR(hydroc::getDataDir());

    // auto h5fname = (DATADIR / "sphere" / "hydroData" /"sphere.h5").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "rm3" / "hydroData" / "rm3.h5").generic_string();
    auto b1Meshfname = (DATADIR / "rm3" / "geometry" / "float_cog.obj").generic_string();
    auto b2Meshfname = (DATADIR / "rm3" / "geometry" / "plate_cog.obj").generic_string();

    // std::cout << "Attempting to open mesh file: " <<
    // std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/meshFiles/float.obj").c_str()) << std::endl;
    double density = 0.0;
    bool evaluate_mass = false;
    bool create_visu_mesh = true;
    bool detect_collision = true;

    auto body1 =
        chrono_types::make_shared<chrono::ChBodyEasyMesh>(b1Meshfname,  // file name
                                                          density, evaluate_mass, create_visu_mesh, detect_collision);
    auto body2 =
        chrono_types::make_shared<chrono::ChBodyEasyMesh>(b2Meshfname,  // file name
                                                          density, evaluate_mass, create_visu_mesh, detect_collision);

    std::vector<H5FileInfo> infos;
    infos.emplace_back(h5fname, "body1");
    infos.emplace_back(h5fname, "body2");  // reopen the file

    std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;

    /// TODO Check if local vector is really copied into constructor of ChLoadAddedMass
    /// else it could be a memory fault

    const size_t nBodies = 2;
    std::vector<std::shared_ptr<ChLoadable>> loadables;
    loadables.push_back(body1);
    loadables.push_back(body2);

    my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(infos, loadables);

    std::cout << "End" << std::endl;
    return 0;
}