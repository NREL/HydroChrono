#include <cstdio>
#include <filesystem>

#include <chrono/solver/ChSolverPMINRES.h>
#include <chrono/solver/ChIterativeSolverLS.h>
#include <chrono/timestepper/ChTimestepper.h>
#include <chrono/timestepper/ChTimestepperHHT.h>

#include <chrono/physics/ChForce.h>
#include <chrono/physics/ChLoadContainer.h>
#include <chrono/physics/ChLoadsBody.h>
#include <chrono/physics/ChLoad.h>
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/physics/ChSystemSMC.h>
#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChBodyEasy.h>

#include <chrono/fea/ChMeshFileLoader.h>


#include <hydroc/h5fileinfo.h>

using namespace chrono;
using namespace chrono::fea;

enum class WaveMode {
    /// @brief No waves
    noWaveCIC = 0,
    /// @brief Regular waves
    regular   = 1,
	irregular = 2
}; 

// =============================================================================
struct HydroInputs {

	WaveMode mode;
    HydroInputs();
    void UpdateNumTimesteps();
    void UpdateRampTimesteps();
    void CreateSpectrum();
    void CreateFreeSurfaceElevation();
    std::vector<double> spectrum_frequencies;
    std::vector<double> spectral_densities;
    double ramp_duration;
	double freq_index_des;
    double wave_height;
    double wave_period;
    double simulation_duration;
    double simulation_dt;
    int num_timesteps;
    int ramp_timesteps;
    std::vector<double> ramp;
	double regular_wave_amplitude;
	double regular_wave_omega;
	double wave_omega_delta;
	std::vector<double> excitation_force_mag;
	std::vector<double> excitation_force_phase;
	HydroInputs(HydroInputs& old) = default;
	HydroInputs& operator = (const HydroInputs& rhs) = default;

};

// =============================================================================
class ForceFunc6d;

class TestHydro;

class ComponentFunc : public ChFunction {
public:
	ComponentFunc();
	ComponentFunc(const ComponentFunc& old);
	ComponentFunc(ForceFunc6d* b, int i);
	virtual ComponentFunc* Clone() const override;
	virtual double 	Get_y(double x) const override;
private:
	ForceFunc6d* base; 
	int index;
};

// =============================================================================
// ForceFunc6d organizes the functional (time dependent) forces in each DoF (6 total) for a body
class ForceFunc6d {
public:
	ForceFunc6d();
	ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* all_hydro_forces_user);
	ForceFunc6d(const ForceFunc6d& old);
	double coordinateFunc(int i);
private:
	void SetForce();
	void SetTorque();
	void ApplyForceAndTorqueToBody();
	std::shared_ptr<ChBody> body;
	int b_num;
	ComponentFunc forces[6];
	std::shared_ptr<ComponentFunc> force_ptrs[6];
	std::shared_ptr<ChForce> chrono_force;
	std::shared_ptr<ChForce> chrono_torque;
	TestHydro* all_hydro_forces;
};

class ChLoadAddedMass;

class TestHydro {
public:
	bool printed = false;
	TestHydro();
	TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies, std::string h5_file_name, HydroInputs& users_hydro_inputs);
	TestHydro(const TestHydro& old) = delete;
	TestHydro operator = (const TestHydro& rhs) = delete;
	void WaveSetUp();
	std::vector<double> ComputeForceHydrostatics();
	std::vector<double> ComputeForceRadiationDampingConv(); 
	std::vector<double> ComputeForceExcitationRegularFreq();
	//std::vector<double> ComputeForceRegularWaves();
	double GetRIRFval(int row, int col, int st);
	double coordinateFunc(int b, int i);
	bool convTrapz;
private:
	std::vector<std::shared_ptr<ChBody>> bodies;
	std::vector<H5FileInfo> file_info;
	std::vector<ForceFunc6d> force_per_body;
	double sumVelHistoryAndRIRF;
	HydroInputs hydro_inputs;
	std::vector<double> force_hydrostatic;
	std::vector<double> force_radiation_damping;
	std::vector<double> force_excitation_freq;
	//std::vector<double> force_reg_waves;
	std::vector<double> total_force;
	int num_bodies;
	std::vector<double> equilibrium;
	std::vector<double> cb_minus_cg;
	double rirf_timestep;
	double getVelHistoryVal(int step, int c) const;
	double setVelHistory(double val, int step, int b_num, int index);

	//double freq_index_des;
	//int freq_index_floor;
	//double freq_interp_val;
	std::vector<double> velocity_history; // use helper function to access vel_history elements correctly
	double prev_time;
	std::vector<double> rirf_time_vector; // (should be the same for each body?)
	int offset_rirf;
	std::shared_ptr<ChLoadContainer> my_loadcontainer;
	std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
};
