#include <hydroc/wave_types.h>

// NoWave class definitions:
Eigen::VectorXd NoWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies * 6;
    Eigen::VectorXd f(dof);
    f.Zero();
    return f;
}

// Regular wave class definitions:
RegularWave::RegularWave() {
    num_bodies = 1;
}

RegularWave::RegularWave(unsigned int num_b) {
    num_bodies = num_b;
}

//void RegularWave::Initialize() {
//    // set up regular waves here, call other helper functions as necessary 
//    excitation_force_mag.resize(total_dofs, 0.0);
//    excitation_force_phase.resize(total_dofs, 0.0);
//    force_excitation_freq.resize(total_dofs, 0.0);
//    wave_omega_delta = file_info[0].GetOmegaDelta();
//    freq_index_des   = (regular_wave_omega / wave_omega_delta) - 1;
//    for (int b = 0; b < num_bodies; b++) {
//        for (int rowEx = 0; rowEx < 6; rowEx++) {
//            int body_offset = 6 * b;
//            excitation_force_mag[body_offset + rowEx] =
//                file_info[b].GetExcitationMagInterp(rowEx, 0, freq_index_des);
//            excitation_force_phase[body_offset + rowEx] =
//                file_info[b].GetExcitationPhaseInterp(rowEx, 0, freq_index_des);
//        }
//    }
//}

Eigen::VectorXd RegularWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies * 6;
    Eigen::VectorXd f(dof);
    // initialize the force here:
    // for now set f to all zeros
    f.Zero();
    return f;
}

// Irregular wave class definitions:
IrregularWave::IrregularWave() {
    num_bodies = 1;
}

IrregularWave::IrregularWave(unsigned int num_b) {
    num_bodies = num_b;
}

void IrregularWave::Initialize() {
    // set up regular waves here, call other helper functions as necessary
}

Eigen::VectorXd IrregularWave::GetForceAtTime(double t) {
    unsigned int dof = num_bodies * 6;
    Eigen::VectorXd f(dof);
    // initialize the force here:
    // for now set f to all zeros
    f.Zero();
    return f;
}