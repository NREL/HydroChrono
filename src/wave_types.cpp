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
        }  // hwy is excitation_force_phase nullptr?
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
    // call resample, call create spectrum and freesurface functions
    // ex_irf
    std::vector<Eigen::MatrixXd> ex_irf_old(num_bodies);
    std::vector<Eigen::VectorXd> ex_irf_time_old(num_bodies);
    // excitation irf same for each body?
    for (int b = 0; b < num_bodies; b++) {
        ex_irf_old[b]      = GetExcitationIRF(b);
        ex_irf_time_old[b] = info[b].excitation_irf_time;
    }

    // resample excitation IRF time series
    // h5 file irf has different timestep, want to resample with interpolation (cubic spline?) once at start no
    // interpolation in convolution integral part
    // different one for each body it's 6x1x1000 so maybe switch to 2d reading
    std::vector<Eigen::MatrixXd> ex_irf_resampled(num_bodies);
    std::vector<Eigen::VectorXd> ex_irf_time_resampled(num_bodies);
    for (int b = 0; b < num_bodies; b++) {  // function call (time_in, vals_in, &time_out, &vals_out)
        ex_irf_time_resampled[b] = ResampleTime(ex_irf_time_old[b], simulation_dt);
        ResampleVals(ex_irf_time_old[b], ex_irf_old[b], ex_irf_time_resampled[b]);
    }
    // vvv this works
    //Eigen::MatrixXd ex_irf_resampled;
    //Eigen::VectorXd ex_irf_time_resampled = ResampleTime(ex_irf_time_old[0], simulation_dt);
     //   ResampleVals(ex_irf_time_old[0], ex_irf_old[0], ex_irf_time_resampled);


    // Eigen::VectorXd ex_irf_resampled = file_info[0].GetExcitationIRFResampled();
    // WriteContainerToFile(ex_irf_resampled, "ex_irf_resampled.txt");
}

#include <unsupported/Eigen/Splines>

Eigen::MatrixXd IrregularWave::ResampleVals(const Eigen::VectorXd& t_old,
                                            Eigen::MatrixXd& vals_old,
                                            const Eigen::VectorXd& t_new) {
    assert(vals_old.rows() == 6);
    
    Eigen::MatrixXd vals_new(6, t_new.size());
    // we need to ensure the old times used start at 0, since the new times will
    double dt_old                 = t_old[1] - t_old[0];
    Eigen::VectorXd t_old_shifted = Eigen::VectorXd::LinSpaced(t_old.size(), 0, (t_old.size() - 1) * dt_old);

    // interpolate the irf over each dof separately, 1 row at a time
    for (int dof = 0; dof < 6; dof++) {
        Eigen::Spline<double, 1> spline =
            Eigen::SplineFitting<Eigen::Spline<double, 1>>::Interpolate(vals_old.row(dof), 3, t_old_shifted);
        for (int i = 0; i < t_new.size(); i++) {
            vals_new(dof, i) = spline(t_new[i])[0];
        }
    }

    std::ofstream out("resample.txt");
    // print for testing if it gets this far:
    for (int i = 0; i < t_new.size(); i++) {
        out << t_new[i];
        for (int dof = 0; dof < 6; dof++) {
            out << " " << vals_new(dof, i);
        }
        out << std::endl;
    }
    out.close();
    out.open("compare.txt");
    for (int i = 0; i < t_old.size(); i++) {
        out << t_old[i];
        for (int dof = 0; dof < 6; dof++) {
            out << " " << vals_old(dof, i);
        }
        out << std::endl;
    }
    out.close();

    return vals_new;
}

Eigen::VectorXd IrregularWave::ResampleTime(const Eigen::VectorXd& t_old, const double dt_new) {
    double dt_old  = t_old[1] - t_old[0];
    int size_new   = static_cast<int>(ceil(t_old.size() * dt_old / dt_new));
    double t_final = (t_old.size() - 1) * dt_old;

    Eigen::VectorXd t_new = Eigen::VectorXd::LinSpaced(size_new, 0, t_final);

    // print for testing
    int N = t_old.size() - 1;
    std::cout << "T_old = {t_i | i = 0 .. " << N << ", t_0 = " << t_old[0] << ", t_" << N << " = " << t_old[N] << "}\n"
              << "dt_old = " << dt_old << "\nT_new = {t_i | i = 0 .. " << N << ", t_0 = " << t_new[0] << ", t_" << N
              << " = " << t_new[N] << "}\n"
              << "dt_new = " << dt_new << std::endl;

        return t_new;
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
    // see ComputeForceExcitation and convolution functions

    /*******************************************************************************
     * TestHydro::ComputeForceExcitation()
     * computes the 6N dimensional excitation force
     *******************************************************************************/
    // TODO delete this
    // std::vector<double> TestHydro::ComputeForceExcitation() {
    //    double time = bodies[0]->GetChTime();
    //
    //    int total_dofs = 6 * num_bodies;
    //    force_excitation.resize(total_dofs, 0.0);
    //
    //    for (int body = 0; body < num_bodies; body++) {
    //        // Loop through the DOFs
    //        for (int dof = 0; dof < 6; ++dof) {
    //            // Compute the convolution for the current DOF
    //            double force_excitation_dof =
    //                ExcitationConvolution(body, dof, time, hydro_inputs.eta, t_irf, hydro_inputs.simulation_dt);
    //            int force_excitation_index               = body * 6 + dof;
    //            force_excitation[force_excitation_index] = force_excitation_dof;
    //        }
    //    }
    //    return force_excitation;
    //}

    return f;
}

/*******************************************************************************
 * IrregularWave::GetExcitationIRF()
 * returns the std::vector of excitation_irf_matrix from h5 file
 *******************************************************************************/
Eigen::MatrixXd IrregularWave::GetExcitationIRF(int b) const {
    return info[b].excitation_irf_matrix;
}