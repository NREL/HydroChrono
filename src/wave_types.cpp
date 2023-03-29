#include <hydroc/wave_types.h>

// NoWave class definitions:
Eigen::VectorXd NoWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies * 6;
    Eigen::VectorXd f(dof);
    for (int i = 0; i < dof; i++) {
        f[i] = 0.0;
    }
    return f;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Regular wave class definitions:
RegularWave::RegularWave() {
    num_bodies = 1;
}

RegularWave::RegularWave(unsigned int num_b) {
    num_bodies = num_b;
}

void RegularWave::Initialize() {
    // set up regular waves here, call other helper functions as necessary
    int total_dofs = 6 * num_bodies;
    excitation_force_mag.resize(total_dofs);
    excitation_force_phase.resize(total_dofs);
    force.resize(total_dofs);

    double wave_omega_delta = GetOmegaDelta();
    double freq_index_des   = (regular_wave_omega / wave_omega_delta) - 1;
    for (int b = 0; b < num_bodies; b++) {
        for (int rowEx = 0; rowEx < 6; rowEx++) {
            int body_offset = 6 * b;
            // why are these always 0? vvv TODO check/change this
            excitation_force_mag[body_offset + rowEx]   = GetExcitationMagInterp(b, rowEx, 0, freq_index_des);
            excitation_force_phase[body_offset + rowEx] = GetExcitationPhaseInterp(b, rowEx, 0, freq_index_des);
        }
    }
}

void RegularWave::AddH5Data(std::vector<HydroData::RegularWaveInfo>& reg_h5_data) {
    info = reg_h5_data;
}

// should return a 6N long vector
Eigen::VectorXd RegularWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies * 6;
    Eigen::VectorXd f(dof);
    // initialize the force here:
    for (int b = 0; b < num_bodies; b++) {
        int body_offset = 6 * b;
        for (int rowEx = 0; rowEx < 6; rowEx++) {
            f[body_offset + rowEx] = excitation_force_mag[body_offset + rowEx] * regular_wave_amplitude *
                                     cos(regular_wave_omega * t + excitation_force_phase[rowEx]);
        } // hwy is excitation_force_phase nullptr?
    }
    return f;
}

// put more reg wave forces here:
// helper GetOmegaDelta()
/*******************************************************************************
 * RegularWave::GetOmegaDelta()
 * returns omega step size
 *******************************************************************************/
double RegularWave::GetOmegaDelta() const {
    double omega_max = info[0].freq_list[info[0].freq_list.size() - 1];
    double num_freqs = info[0].freq_list.size();
    return omega_max / num_freqs;
}

/*******************************************************************************
 * RegularWave::GetExcitationMagInterp()
 * returns excitation magnitudes for body b, row i, column j, frequency ix k
 *******************************************************************************/
double RegularWave::GetExcitationMagInterp(int b, int i, int j, double freq_index_des) const {
    double freq_interp_val    = freq_index_des - floor(freq_index_des);
    double excitationMagFloor = info[b].excitation_mag_matrix(i, j, (int)floor(freq_index_des));
    double excitationMagCeil  = info[b].excitation_mag_matrix(i, j, (int)floor(freq_index_des) + 1);
    double excitationMag      = (freq_interp_val * (excitationMagCeil - excitationMagFloor)) + excitationMagFloor;

    return excitationMag;
}

/*******************************************************************************
 * RegularWave::GetExcitationPhaseInterp()
 * returns excitation phases for row i, column j, frequency ix k
 *******************************************************************************/
double RegularWave::GetExcitationPhaseInterp(int b, int i, int j, double freq_index_des) const {
    double freq_interp_val      = freq_index_des - floor(freq_index_des);  // look into c++ modf TODO
    double excitationPhaseFloor = info[b].excitation_phase_matrix(
        i, j, (int)floor(freq_index_des));  // TODO check if freq_index_des is >0, if so just cast instead of floor
    double excitationPhaseCeil = info[b].excitation_phase_matrix(i, j, (int)floor(freq_index_des) + 1);
    double excitationPhase = (freq_interp_val * (excitationPhaseCeil - excitationPhaseFloor)) + excitationPhaseFloor;

    return excitationPhase;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Irregular wave class definitions:
IrregularWave::IrregularWave() {
    num_bodies = 1;
}

IrregularWave::IrregularWave(unsigned int num_b) {
    num_bodies = num_b;
}

void IrregularWave::Initialize() {
    // set up irregular waves here, call other helper functions as necessary
}

void IrregularWave::AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data) {
    info = irreg_h5_data;
}
Eigen::VectorXd IrregularWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies * 6;
    Eigen::VectorXd f(dof);
    // initialize the force here:
    // for now set f to all zeros
    for (int i = 0; i < dof; i++) {
        f[i] = 0.0;
    }
    return f;
}