#include <hydroc/gui/guihelper.h>
#include <hydroc/helper.h>
#include <hydroc/hydro_system.h>
#include <hydroc/logging.h>

#ifdef HYDROCHRONO_HAVE_MOORDYN
#include <hydroc/config/moordyn_config.h>
#endif

#include <chrono/core/ChRealtimeStep.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChSystemNSC.h>

#include <chrono>
#include <iomanip>
#include <vector>

using namespace chrono;

int main(int argc, char* argv[]) {
    std::cout << "Chrono version: " << CHRONO_VERSION << "\n\n";

    bool profilingOn     = true;
    bool saveDataOn      = true;
    bool plotOn          = true;
    bool visualizationOn = true;
    std::string data_dir;
    if (!hydroc::GetCLIArguments(argc, argv, "RM3 mooring demo", saveDataOn, profilingOn, plotOn, visualizationOn,
                                 data_dir))
        return 1;
    if (!hydroc::SetInitialEnvironment(data_dir)) return 1;

    std::filesystem::path DATADIR(hydroc::getDataDir());
    std::cout << "Data directory: " << DATADIR << std::endl;

    auto body1_meshfame =
        (DATADIR / "demos" / "rm3" / "geometry" / "float_cog.obj").lexically_normal().generic_string();
    auto body2_meshfame =
        (DATADIR / "demos" / "rm3" / "geometry" / "plate_cog.obj").lexically_normal().generic_string();
    auto h5fname = (DATADIR / "demos" / "rm3" / "hydroData" / "rm3.h5").lexically_normal().generic_string();

    std::cout << "  body1 mesh: " << body1_meshfame << std::endl;
    std::cout << "  body2 mesh: " << body2_meshfame << std::endl;
    std::cout << "  h5 file:    " << h5fname << std::endl;

    // System and solver
    ChSystemNSC system;
    system.SetGravitationalAcceleration(ChVector3d(0.0, 0.0, -9.81));
    double timestep = 0.01;
    system.SetTimestepperType(ChTimestepper::Type::HHT);
    system.SetSolverType(ChSolver::Type::SPARSE_LU);
    ChRealtimeStepTimer realtime_timer;
    double simulationDuration = 40.0;

    std::shared_ptr<hydroc::gui::UI> pui = hydroc::gui::CreateUI(visualizationOn);
    hydroc::gui::UI& ui = *pui.get();

    std::vector<double> time_vector;
    std::vector<double> float_heave_position;
    std::vector<double> float_drift_position;
    std::vector<double> plate_heave_position;

    // Float body
    std::cout << "Loading float mesh: " << body1_meshfame << std::endl;
    std::shared_ptr<ChBody> float_body1 = chrono_types::make_shared<ChBodyEasyMesh>(
        body1_meshfame, 0, false, true, false);

    system.Add(float_body1);
    float_body1->SetName("body1");
    float_body1->SetPos(ChVector3d(0, 0, -0.72));
    float_body1->SetMass(725834);
    float_body1->SetInertiaXX(ChVector3d(20907301.0, 21306090.66, 37085481.11));

    auto red = chrono_types::make_shared<ChVisualMaterial>();
    red->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.1f));
    float_body1->GetVisualShape(0)->SetMaterial(0, red);

    // Plate body
    std::cout << "Loading plate mesh: " << body2_meshfame << std::endl;
    std::shared_ptr<ChBody> plate_body2 = chrono_types::make_shared<ChBodyEasyMesh>(
        body2_meshfame, 0, false, true, false);

    system.Add(plate_body2);
    plate_body2->SetName("body2");
    plate_body2->SetPos(ChVector3d(0, 0, -21.29));
    plate_body2->SetMass(886691);
    plate_body2->SetInertiaXX(ChVector3d(94419614.57, 94407091.24, 28542224.82));

    auto blue = chrono_types::make_shared<ChVisualMaterial>();
    blue->SetDiffuseColor(ChColor(0.3f, 0.1f, 0.6f));
    plate_body2->GetVisualShape(0)->SetMaterial(0, blue);

    // Prismatic joint (vertical heave)
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(float_body1, plate_body2, false, ChFramed(ChVector3d(0, 0, -0.72)),
                          ChFramed(ChVector3d(0, 0, -21.29)));
    system.AddLink(prismatic);

    // PTO damper
    auto prismatic_pto = chrono_types::make_shared<ChLinkTSDA>();
    prismatic_pto->Initialize(float_body1, plate_body2, false, ChVector3d(0, 0, -0.72), ChVector3d(0, 0, -21.29));
    prismatic_pto->SetDampingCoefficient(1200000.0);
    system.AddLink(prismatic_pto);

    // Hydrodynamic forces
    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(float_body1);
    bodies.push_back(plate_body2);

    auto my_hydro_inputs                     = std::make_shared<RegularWave>();
    my_hydro_inputs->regular_wave_amplitude_ = 1.0;
    my_hydro_inputs->regular_wave_omega_     = 2.10;

    TestHydro hydro_forces(bodies, h5fname);

#ifdef HYDROCHRONO_HAVE_MOORDYN
    {
        auto moordyn_input =
            (DATADIR / "yaml" / "rm3" / "mooring" / "lines_rm3.txt").lexically_normal().generic_string();
        std::cout << "  MoorDyn input: " << moordyn_input << std::endl;

        MoorDynConfig md_cfg;
        md_cfg.enabled = true;
        md_cfg.input_file = moordyn_input;
        md_cfg.coupled_body_indices = {1};  // body2 (plate) = index 1
        hydro_forces.SetMoorDynConfig(md_cfg);
    }
#endif

    hydro_forces.AddWaves(my_hydro_inputs);

    auto start = std::chrono::high_resolution_clock::now();

    ui.Init(&system, "RM3 - Mooring Demo");
    ui.SetCamera(0, -50, -10, 0, 0, -10);
    ui.SetWaveModel(my_hydro_inputs);

    while (system.GetChTime() <= simulationDuration) {
        if (ui.IsRunning(timestep) == false) break;

        if (ui.simulationStarted) {
            system.DoStepDynamics(timestep);

            time_vector.push_back(system.GetChTime());
            float_heave_position.push_back(float_body1->GetPos().z());
            float_drift_position.push_back(float_body1->GetPos().x());
            plate_heave_position.push_back(plate_body2->GetPos().z());
        }
    }

    auto end          = std::chrono::high_resolution_clock::now();
    unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::string out_dir = hydroc::getDemoOutDir();
    if (profilingOn || saveDataOn) {
        out_dir = out_dir + "/" + RESULTS_DIR_NAME;
        std::filesystem::create_directory(std::filesystem::path(out_dir));
    }

    if (profilingOn) {
        std::ofstream profilingFile(out_dir + "/mooring_duration.txt");
        profilingFile << duration << " ms\n";
        profilingFile.close();
    }

    if (saveDataOn) {
        std::ofstream outputFile(out_dir + "/mooring.txt");
        outputFile << std::left << std::setw(10) << "Time (s)" << std::right << std::setw(16) << "Float Heave (m)"
                   << std::right << std::setw(16) << "Plate Heave (m)" << std::right << std::setw(16)
                   << "Float Drift (x) (m)" << std::endl;
        for (size_t i = 0; i < time_vector.size(); ++i)
            outputFile << std::left << std::setw(10) << std::setprecision(2) << std::fixed << time_vector[i]
                       << std::right << std::setw(16) << std::setprecision(8) << std::fixed << float_heave_position[i]
                       << std::right << std::setw(16) << std::setprecision(8) << std::fixed << plate_heave_position[i]
                       << std::right << std::setw(16) << std::setprecision(8) << std::fixed << float_drift_position[i]
                       << std::endl;
        outputFile.close();
    }
    return 0;
}
