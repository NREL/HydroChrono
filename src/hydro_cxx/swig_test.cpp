#include "swig_test.h"
#include <iostream>

void test() {
	std::cout << "Hello World!" << std::endl;
}

// =============================================================================
// HydroInputs Class Definitions
// =============================================================================

/*******************************************************************************
* HydroInputs constructor
*******************************************************************************/
HydroInputs::HydroInputs() {
	// TODO: switch depending on wave option (regular, regularCIC, irregular, noWaveCIC) enum?
	mode = NONE;
	regular_wave_amplitude = 0;
	//excitation_force_phase.resize(6,0);
	//excitation_force_mag.resize(6,0);
}

/*******************************************************************************
* HydroInputs constructor
*******************************************************************************/
HydroInputs::HydroInputs(HydroInputs& old) {
	*this = old;
}

/*******************************************************************************
* HydroInputs constructor
*******************************************************************************/
HydroInputs& HydroInputs::operator = (HydroInputs& rhs) {
	freq_index_des = rhs.freq_index_des;
	regular_wave_amplitude = rhs.regular_wave_amplitude;
	regular_wave_omega = rhs.regular_wave_omega;
	wave_omega_delta = rhs.wave_omega_delta;
	excitation_force_mag = rhs.excitation_force_mag;
	excitation_force_phase = rhs.excitation_force_phase;
	mode = rhs.mode;
	return *this;
}
