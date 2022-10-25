#ifndef _SWIG_TEST_
#define _SWIG_TEST_
#include <vector>
void test();
enum WaveMode { NONE, REGULAR }; // eventually add irregular waves mode
// =============================================================================
class HydroInputs {
public:
	WaveMode mode;
	HydroInputs();
	double freq_index_des;
	double regular_wave_amplitude;
	double regular_wave_omega;
	double wave_omega_delta;
	std::vector<double> excitation_force_mag;
	std::vector<double> excitation_force_phase;
	HydroInputs(HydroInputs& old);
	HydroInputs& operator = (HydroInputs& rhs);
private:
};
#endif // !