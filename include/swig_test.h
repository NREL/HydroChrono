#ifndef _SWIG_TEST_
#define _SWIG_TEST_

#ifndef SWIG
#include <vector>
#endif

#include "hydro_forces.h"

void test();
enum WaveMode { NONE, REGULAR }; // eventually add irregular waves mode
// =============================================================================
class HydroInputs {
public:
	WaveMode mode;
	HydroInputs();
	~HydroInputs() = default;
	double freq_index_des;
	double regular_wave_amplitude;
	double regular_wave_omega;
	double wave_omega_delta;
	std::vector<double> excitation_force_mag;
	std::vector<double> excitation_force_phase;
	HydroInputs(const HydroInputs& old);
	HydroInputs& operator = (const HydroInputs& rhs);
	void test2();
private:
};

class ChLoadAddedMass;
class H5FileInfo;
class ForceFunc6d;
class TestHydro {
public:
	bool printed = false;
	TestHydro();
	TestHydro(const std::vector<std::shared_ptr<chrono::ChBody>>& user_bodies, std::string h5_file_name, HydroInputs& users_hydro_inputs);
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
	std::vector<std::shared_ptr<chrono::ChBody>> bodies;
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
	double getVelHistoryAllBodies(int step, int c) const;
	double setVelHistory(double val, int step, int b_num, int index);

	//double freq_index_des;
	//int freq_index_floor;
	//double freq_interp_val;
	std::vector<double> velocity_history; // use helper function to access vel_history elements correctly
	double prev_time;
	std::vector<double> rirf_time_vector; // (should be the same for each body?)
	int offset_rirf;
	std::shared_ptr<chrono::ChLoadContainer> my_loadcontainer;
	std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
};

#endif // !