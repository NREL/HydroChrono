#include <hydroc/chloadaddedmass.h>
#include <hydroc/h5fileinfo.h>
#include <hydroc/helper.h>

#include <chrono/core/ChTypes.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemSMC.h>

#include <cstdlib>
#include <filesystem>  // C++17
#include <iostream>
#include <memory>
#include <vector>

using std::filesystem::absolute;
using std::filesystem::path;

int main(int argc, char* argv[]) {
    if (hydroc::SetInitialEnvironment(argc, argv) != 0) {
        return 1;
    }

    path DATADIR(hydroc::getDataDir());

    // auto h5fname = (DATADIR / "sphere" / "hydroData" /"sphere.h5").lexically_normal().generic_string();
    auto h5fname     = (DATADIR / "rm3" / "hydroData" / "rm3.h5").generic_string();
    auto b1Meshfname = (DATADIR / "rm3" / "geometry" / "float_cog.obj").generic_string();
    auto b2Meshfname = (DATADIR / "rm3" / "geometry" / "plate_cog.obj").generic_string();

    // std::cout << "Attempting to open mesh file: " <<
    // std::filesystem::absolute(GetChronoDataFile("../../HydroChrono/meshFiles/float.obj").c_str()) << std::endl;
    double density        = 0.0;
    bool evaluate_mass    = false;
    bool create_visu_mesh = true;
    bool detect_collision = true;

    auto body1 =
        chrono_types::make_shared<chrono::ChBodyEasyMesh>(b1Meshfname,  // file name
                                                          density, evaluate_mass, create_visu_mesh, detect_collision);
    auto body2 =
        chrono_types::make_shared<chrono::ChBodyEasyMesh>(b2Meshfname,  // file name
                                                          density, evaluate_mass, create_visu_mesh, detect_collision);

    HydroData infos = H5FileInfo(h5fname, 2).ReadH5Data();

    std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;

    /// TODO Check if local vector is really copied into constructor of ChLoadAddedMass
    /// else it could be a memory fault

    const size_t nBodies = 2;
    std::vector<std::shared_ptr<ChLoadable>> loadables;
    loadables.push_back(body1);
    loadables.push_back(body2);

    ChSystemSMC my_system;

    my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(infos.GetBodyInfos(), loadables, &my_system);

    std::cout << "End" << std::endl;
    return 0;
}