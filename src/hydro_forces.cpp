/*********************************************************************
 * @file  hydro_forces.cpp
 *
 * @brief Implementation of TestHydro main class and helper classes
 * ComponentFunc and ForceFunc6d.
 *
 * OVERVIEW:
 * This file implements the core hydrodynamic force computation for multibody
 * marine systems. It connects Chrono physics bodies to hydrodynamic models
 * (hydrostatics, radiation damping, wave excitation) and applies the
 * resulting forces during simulation.
 *
 * MAIN RESPONSIBILITIES:
 * - TestHydro: Orchestrates all hydrodynamic force components for multiple bodies
 * - ForceFunc6d: Wraps Chrono ChForce/ChTorque callbacks for each body
 * - ComponentFunc: Provides per-DOF force values to Chrono's force system
 *
 * INTERACTIONS:
 * - Reads HDF5 hydrodynamic data via H5FileInfo (equilibrium, stiffness, RIRF kernels)
 * - Uses WaveBase hierarchy for wave excitation forces
 * - Applies forces to Chrono bodies through ChForce/ChLoadAddedMass
 * - Maintains velocity history buffers for radiation convolution
 *
 * KEY ASSUMPTIONS:
 * - All bodies share the same H5 file (multibody data in single file)
 * - Bodies are 1-indexed in ForceFunc6d interface (legacy)
 * - 6 DOF per body (surge, sway, heave, roll, pitch, yaw)
 * - Forces computed once per time step, cached via prev_time check
 * - Radiation history stored per-body in velocity_history_ vector
 *
 * KNOWN LIMITATIONS:
 * - Monolithic design: all force models mixed in one class
 * - Tight coupling to Chrono types (hard to test without Chrono)
 * - Body indexing inconsistency (1-indexed in some places, 0-indexed in others)
 * - No per-body enable/disable of radiation or excitation (system-wide only)
 *
 * DEBUG INSTRUMENTATION:
 * To enable debug output, compile with -DHYDROCHRONO_DEBUG.
 * This will print force components, history sizes, and timing info.
 * Baseline outputs to record:
 * - Total force vector per time step (total_force_)
 * - Individual components (hydrostatic, radiation, waves)
 * - Velocity history size and time span
 * - RIRF kernel access patterns
 *********************************************************************/

// TODO minimize include statements, move all to header file hydro_forces.h?
#include "hydroc/hydro_forces.h"
#include <hydroc/chloadaddedmass.h>
#include <hydroc/h5fileinfo.h>
#include <hydroc/wave_types.h>
#include <hydroc/logging.h>
#include "hydro/core/system_state.h"
#include "hydro/radiation/radiation_rirf_processing.h"
#include "hydro/radiation/radiation_rirf_convolution.h"

#include <chrono/physics/ChLoad.h>
#include <unsupported/Eigen/Splines>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>  // std::accumulate
#include <random>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <filesystem>

#ifdef _OPENMP
#include <omp.h>
#endif

const int kDofPerBody  = 6;
const int kDofLinOrRot = 3;

// ------------------------------------------------------------
// SECTION: Utility functions (legacy)
// ------------------------------------------------------------

// Internal helper: Build SystemState from all Chrono bodies (legacy implementation).
// This extracts pose and velocity data from Chrono bodies into the new SystemState structure.
namespace {
hydrochrono::hydro::SystemState BuildSystemStateFromChronoBodies(
    const std::vector<std::shared_ptr<ChBody>>& bodies) {
    hydrochrono::hydro::SystemState system_state;
    system_state.bodies.reserve(bodies.size());
    
    for (const auto& body : bodies) {
        hydrochrono::hydro::BodyState body_state;
        
        // Extract position and orientation
        const chrono::ChVector3d position_world = body->GetPos();
        const chrono::ChVector3d rotation_rpy    = body->GetRot().GetCardanAnglesXYZ();
        
        body_state.position[0] = position_world.x();
        body_state.position[1] = position_world.y();
        body_state.position[2] = position_world.z();
        
        body_state.orientation_rpy[0] = rotation_rpy.x();
        body_state.orientation_rpy[1] = rotation_rpy.y();
        body_state.orientation_rpy[2] = rotation_rpy.z();
        
        // Extract velocities
        const auto linear_velocity_world  = body->GetPosDt();
        const auto angular_velocity_world = body->GetAngVelParent();
        
        body_state.linear_velocity[0] = linear_velocity_world.x();
        body_state.linear_velocity[1] = linear_velocity_world.y();
        body_state.linear_velocity[2] = linear_velocity_world.z();
        
        body_state.angular_velocity[0] = angular_velocity_world.x();
        body_state.angular_velocity[1] = angular_velocity_world.y();
        body_state.angular_velocity[2] = angular_velocity_world.z();
        
        system_state.bodies.push_back(body_state);
    }
    
    return system_state;
}
}  // namespace

/**
 * @brief Generates a vector of evenly spaced numbers over a specified range.
 *
 * This function returns a vector of `num_points` numbers evenly spaced from
 * `start` to `end`. The function utilizes a single loop for this computation,
 * making it efficient for generating large vectors.
 *
 * @param start - The start value of the sequence.
 * @param end - The end value of the sequence.
 * @param num_points - The number of evenly spaced samples to generate.
 * @return std::vector<double> - Vector of evenly spaced numbers.
 * @exception None
 */
std::vector<double> Linspace(double start, double end, int num_points) {
    std::vector<double> result(num_points);
    double step = (end - start) / (num_points - 1);

    for (int i = 0; i < num_points; ++i) {
        result[i] = start + i * step;
    }

    return result;
}

// ------------------------------------------------------------
// SECTION: Chrono coupling (connects hydrodynamics to ChBody)
// ------------------------------------------------------------
// These classes bridge Chrono's force callback system to TestHydro.
// TODO: Extract to separate module in future refactor.

// TODO reorder ComponentFunc implementation functions to match the header order of functions
ComponentFunc::ComponentFunc() {
    base_  = NULL;
    index_ = kDofPerBody;
}

ComponentFunc::ComponentFunc(ForceFunc6d* b, int i) : base_(b), index_(i) {}

ComponentFunc* ComponentFunc::Clone() const {
    return new ComponentFunc(*this);
}

ComponentFunc::ComponentFunc(const ComponentFunc& old) {
    base_  = old.base_;
    index_ = old.index_;
}

double ComponentFunc::GetVal(double x) const {
    if (base_ == NULL) {
        std::cout << "base == Null!" << std::endl;
        return 0;
    }
    return base_->CoordinateFunc(index_);
}

ForceFunc6d::ForceFunc6d() : forces_{{this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5}} {
    for (unsigned i = 0; i < 6; i++) {
        force_ptrs_[i] = std::shared_ptr<ComponentFunc>(forces_ + i, [](ComponentFunc*) {});
        // sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is
        // automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
        // them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the
        // default deletion logic to do nothing
        // Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
    }
    chrono_force_  = chrono_types::make_shared<ChForce>();
    chrono_torque_ = chrono_types::make_shared<ChForce>();
    chrono_force_->SetAlign(ChForce::AlignmentFrame::WORLD_DIR);
    chrono_torque_->SetAlign(ChForce::AlignmentFrame::WORLD_DIR);
    chrono_force_->SetName("hydroforce");
    chrono_torque_->SetName("hydrotorque");
}

ForceFunc6d::ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* user_all_forces) : ForceFunc6d() {
    body_             = object;
    std::string temp  = body_->GetName();        // remove "body" from "bodyN", convert N to int, get body num
    b_num_            = stoi(temp.erase(0, 4));  // 1 indexed TODO: fix b_num starting here to be 0 indexed
    all_hydro_forces_ = user_all_forces;         // TODO switch to smart pointers? does this use = ?
    if (all_hydro_forces_ == NULL) {
        std::cout << "all hydro forces null " << std::endl;
    }
    SetForce();
    SetTorque();
    ApplyForceAndTorqueToBody();
}

ForceFunc6d::ForceFunc6d(const ForceFunc6d& old)
    : forces_{{this, 0}, {this, 1}, {this, 2}, {this, 3}, {this, 4}, {this, 5}} {
    for (unsigned i = 0; i < 6; i++) {
        force_ptrs_[i] = std::shared_ptr<ComponentFunc>(forces_ + i, [](ComponentFunc*) {});
        // sets force_ptrs[i] to point to forces[i] but since forces is on the stack, it is faster and it is
        // automatically deallocated...shared pointers typically manage heap pointers, and will try deleting
        // them as soon as done. Doesn't work on stack array (can't delete stack arrays), we overload the
        // default deletion logic to do nothing
        // Also! don't need to worry about deleting this later, because stack arrays are always deleted automatically
    }
    chrono_force_     = old.chrono_force_;
    chrono_torque_    = old.chrono_torque_;
    body_             = old.body_;
    b_num_            = old.b_num_;
    all_hydro_forces_ = old.all_hydro_forces_;
    SetForce();
    SetTorque();
}

double ForceFunc6d::CoordinateFunc(int i) {
    // b_num is 1 indexed?
    if (i >= kDofPerBody || i < 0) {
        std::cout << "wrong index force func 6d" << std::endl;
        return 0;
    }
    return all_hydro_forces_->CoordinateFuncForBody(
        b_num_, i);  // b_num is 1 indexed here!!!!! TODO: change all b_num to be 0 indexed everywhere
}

void ForceFunc6d::SetForce() {
    if (chrono_force_ == NULL || body_ == NULL) {
        std::cout << "set force null issue" << std::endl;
    }
    chrono_force_->SetF_x(force_ptrs_[0]);
    chrono_force_->SetF_y(force_ptrs_[1]);
    chrono_force_->SetF_z(force_ptrs_[2]);
}

void ForceFunc6d::SetTorque() {
    if (chrono_torque_ == NULL || body_ == NULL) {
        std::cout << "set torque null issue" << std::endl;
    }
    chrono_torque_->SetF_x(force_ptrs_[3]);
    chrono_torque_->SetF_y(force_ptrs_[4]);
    chrono_torque_->SetF_z(force_ptrs_[5]);
    chrono_torque_->SetMode(ChForce::ForceType::TORQUE);
}

void ForceFunc6d::ApplyForceAndTorqueToBody() {
    body_->AddForce(chrono_force_);
    body_->AddForce(chrono_torque_);
}

// ------------------------------------------------------------
// SECTION: TestHydro class (main hydrodynamics orchestrator)
// ------------------------------------------------------------
// TODO: This class will be refactored into modular force models.
// Current structure mixes all physics components together.

// Explicit destructor definition (needed for unique_ptr to incomplete type)
TestHydro::~TestHydro() = default;

TestHydro::TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
                     std::string h5_file_name,
                     std::shared_ptr<WaveBase> waves)
    : bodies_(user_bodies),
      num_bodies_(bodies_.size()),
      file_info_(H5FileInfo(h5_file_name, num_bodies_).ReadH5Data()) {
    prev_time = -1;

    // Set up time vector
    rirf_time_vector = file_info_.GetRIRFTimeVector();
    // width array
    rirf_width_vector.resize(rirf_time_vector.size());
    for (int ii = 0; ii < rirf_width_vector.size(); ii++) {
        rirf_width_vector[ii] = 0.0;
        if (ii < rirf_time_vector.size() - 1) {
            rirf_width_vector[ii] += 0.5 * abs(rirf_time_vector[ii + 1] - rirf_time_vector[ii]);
        }
        if (ii > 0) {
            rirf_width_vector[ii] += 0.5 * abs(rirf_time_vector[ii] - rirf_time_vector[ii - 1]);
        }
    }

    // Total degrees of freedom (multibody: 6 DOF per body)
    int total_dofs = kDofPerBody * num_bodies_;

    // Initialize radiation convolution module
    const int rirf_steps = file_info_.GetRIRFDims(2);
    radiation_convolution_ = std::make_unique<hydrochrono::hydro::RadiationRirfConvolution>(
        num_bodies_, rirf_steps, rirf_time_vector, rirf_width_vector);

    // Initialize vectors (legacy: kept for compatibility, but radiation_convolution_ owns the history)
    time_history_.clear();
    velocity_history_.clear();
    for (int b = 0; b < num_bodies_; ++b) {
        velocity_history_.push_back(std::vector<std::vector<double>>(0));
    }
    force_hydrostatic_.assign(total_dofs, 0.0);
    force_radiation_damping_.assign(total_dofs, 0.0);
    total_force_.assign(total_dofs, 0.0);
    equilibrium_.assign(total_dofs, 0.0);
    cb_minus_cg_.assign(kDofLinOrRot * num_bodies_, 0.0);

    // Compute equilibrium and cb_minus_cg_ (multibody loop)
    for (int b = 0; b < num_bodies_; ++b) {
        for (int i = 0; i < kDofLinOrRot; ++i) {
            unsigned eq_idx = i + kDofPerBody * b;
            unsigned c_idx  = i + kDofLinOrRot * b;

            equilibrium_[eq_idx] = file_info_.GetCGVector(b)[i];
            cb_minus_cg_[c_idx]  = file_info_.GetCBVector(b)[i] - file_info_.GetCGVector(b)[i];
        }
    }

    // Create Chrono force callbacks for each body (multibody setup)
    for (int b = 0; b < num_bodies_; ++b) {
        force_per_body_.emplace_back(bodies_[b], this);
    }

    // Handle added mass info (applied via Chrono load system)
    my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();

    std::vector<std::shared_ptr<ChLoadable>> loadables(bodies_.size());
    for (int i = 0; i < static_cast<int>(bodies_.size()); ++i) {
        loadables[i] = bodies_[i];
    }

    my_loadbodyinertia =
        chrono_types::make_shared<ChLoadAddedMass>(file_info_.GetBodyInfos(), loadables, bodies_[0]->GetSystem());

    bodies_[0]->GetSystem()->Add(my_loadcontainer);
    my_loadcontainer->Add(my_loadbodyinertia);

    // Set up hydro inputs
    user_waves_ = waves;
    AddWaves(user_waves_);

    // If irregular waves are active, publish spectrum and free-surface inputs into HDF5 if an exporter is present.
    // The exporter is managed by the runner; here we only expose getters via the waves object.
}

void TestHydro::AddWaves(std::shared_ptr<WaveBase> waves) {
    user_waves_ = waves;

    switch (user_waves_->GetWaveMode()) {
        case WaveMode::regular: {
            auto reg = std::static_pointer_cast<RegularWave>(user_waves_);
            reg->AddH5Data(file_info_.GetRegularWaveInfos(), file_info_.GetSimulationInfo());
            break;
        }
        case WaveMode::irregular: {
            auto irreg = std::static_pointer_cast<IrregularWaves>(user_waves_);
            irreg->AddH5Data(file_info_.GetIrregularWaveInfos(), file_info_.GetSimulationInfo());
            break;
        }
    }

    user_waves_->Initialize();
}

// ------------------------------------------------------------
// SECTION: Hydrostatics (legacy)
// ------------------------------------------------------------
// Computes restoring forces from displacement and buoyancy.
// TODO: Extract to HydrostaticsModel class in future refactor.

std::vector<double> TestHydro::ComputeForceHydrostatics() {
    auto __t0 = std::chrono::steady_clock::now();
    assert(num_bodies_ > 0);

    const double rho = file_info_.GetRhoVal();
    const auto gravitational_acceleration = bodies_[0]->GetSystem()->GetGravitationalAcceleration();  // same system for all bodies
    const double rho_times_g = rho * gravitational_acceleration.Length();

#ifdef HYDROCHRONO_DEBUG
    // Debug: log hydrostatics computation start
    std::cout << "[HYDRO_DEBUG] ComputeForceHydrostatics: num_bodies=" << num_bodies_ << ", rho=" << rho << std::endl;
#endif

    // Build SystemState from all Chrono bodies (legacy implementation)
    const auto system_state = BuildSystemStateFromChronoBodies(bodies_);

    // Multibody loop: compute hydrostatic forces for each body
    for (int b = 0; b < num_bodies_; ++b) {
        const auto& body_state = system_state.bodies[b];
        const auto& body = bodies_[b];

        const int body_offset = kDofPerBody * b;
        double* const body_force_hydrostatic = &force_hydrostatic_[body_offset];
        const double* const body_equilibrium = &equilibrium_[body_offset];

        // Current pose (from SystemState)
        const Eigen::Vector3d& position_world = body_state.position;
        const Eigen::Vector3d& rotation_rpy    = body_state.orientation_rpy;

        // 6-DOF displacement from equilibrium (translation xyz, rotation rpy)
        Eigen::Matrix<double, kDofPerBody, 1> displacement_from_equilibrium;
        displacement_from_equilibrium[0] = position_world[0] - body_equilibrium[0];
        displacement_from_equilibrium[1] = position_world[1] - body_equilibrium[1];
        displacement_from_equilibrium[2] = position_world[2] - body_equilibrium[2];
        displacement_from_equilibrium[3] = rotation_rpy[0]   - body_equilibrium[3];
        displacement_from_equilibrium[4] = rotation_rpy[1]   - body_equilibrium[4];
        displacement_from_equilibrium[5] = rotation_rpy[2]   - body_equilibrium[5];

        // Linear hydrostatic restoring force/torque
        const Eigen::MatrixXd restoring_stiffness_matrix = file_info_.GetLinMatrix(b);
        const Eigen::Matrix<double, kDofPerBody, 1> restoring_force_torque =
            -rho_times_g * (restoring_stiffness_matrix * displacement_from_equilibrium);
        for (int i = 0; i < kDofPerBody; ++i) {
            body_force_hydrostatic[i] += restoring_force_torque[i];
        }

        // Buoyancy force at equilibrium: F = rho * (-g) * displaced_volume
        const double displaced_volume = file_info_.GetDispVolVal(b);
        const chrono::ChVector3d buoyancy_force = rho * (-gravitational_acceleration) * displaced_volume;
        body_force_hydrostatic[0] += buoyancy_force.x();
        body_force_hydrostatic[1] += buoyancy_force.y();
        body_force_hydrostatic[2] += buoyancy_force.z();

        // Buoyancy-induced moment about CG: (r_CB - r_CG) x buoyancy
        const int rotation_offset = kDofLinOrRot * b;
        const chrono::ChVector3d cg_to_cb(
            cb_minus_cg_[rotation_offset + 0],
            cb_minus_cg_[rotation_offset + 1],
            cb_minus_cg_[rotation_offset + 2]
        );
        const chrono::ChVector3d buoyancy_torque = cg_to_cb % buoyancy_force;
        body_force_hydrostatic[3] += buoyancy_torque.x();
        body_force_hydrostatic[4] += buoyancy_torque.y();
        body_force_hydrostatic[5] += buoyancy_torque.z();
    }

#ifdef HYDROCHRONO_DEBUG
    // Debug: log hydrostatics result summary
    std::cout << "[HYDRO_DEBUG] Hydrostatics complete: max_force=" 
              << *std::max_element(force_hydrostatic_.begin(), force_hydrostatic_.end()) << std::endl;
#endif

    profile_stats_.hydrostatics_seconds += std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - __t0).count();
    profile_stats_.hydrostatics_calls++;
    return force_hydrostatic_;
}

// ------------------------------------------------------------
// SECTION: Radiation convolution + history buffers
// ------------------------------------------------------------
// Implements velocity-history convolution for radiation damping.
// Maintains time-history buffers per body for multibody systems.
// TODO: Extract to RadiationModel class with dedicated HistoryBuffer.

// Legacy helper functions removed - now in RadiationRirfConvolution class

// Preprocess the radiation kernel K(t) per body for TaperedDirect mode.
void TestHydro::EnsureProcessedRIRF() {
    if (rirf_processed_ready_) {
        return;
    }

    const int steps = file_info_.GetRIRFDims(2);
    
    // Convert TaperedDirectOptions to hydrochrono::hydro::TaperedDirectOptions
    hydrochrono::hydro::TaperedDirectOptions opts;
    opts.smoothing = tapered_opts_.smoothing;
    opts.window_length = tapered_opts_.window_length;
    opts.rirf_end_time = tapered_opts_.rirf_end_time;
    opts.taper_start_percent = tapered_opts_.taper_start_percent;
    opts.taper_end_percent = tapered_opts_.taper_end_percent;
    opts.taper_final_amplitude = tapered_opts_.taper_final_amplitude;
    opts.export_plot_csv = tapered_opts_.export_plot_csv;

    // Create lambda to get RIRF values from file_info_
    auto get_rirf_val = [this](int body, int row_dof, int col, int step) -> double {
        return file_info_.GetRIRFVal(body, row_dof, col, step);
    };

    // Process RIRF kernels using the dedicated module
    rirf_processed_ = hydrochrono::hydro::ProcessRirfKernels(
        num_bodies_, steps, rirf_time_vector, get_rirf_val, opts, diagnostics_output_dir_);

    rirf_processed_ready_ = true;
}

// Main radiation convolution computation.
// Uses velocity history buffers (multibody: one per body).
std::vector<double> TestHydro::ComputeForceRadiationDampingConv() {
    auto __t0 = std::chrono::steady_clock::now();
    const int rirf_steps = file_info_.GetRIRFDims(2);
    const int total_dofs = kDofPerBody * num_bodies_;

    assert(total_dofs > 0 && rirf_steps > 0);

    // If using TaperedDirect, ensure processed kernel is ready
    if (convolution_mode_ == RadiationConvolutionMode::TaperedDirect) {
        EnsureProcessedRIRF();
    }

    // Current time
    const double simulation_time = bodies_[0]->GetChTime();

#ifdef HYDROCHRONO_DEBUG
    std::cout << "[HYDRO_DEBUG] Radiation conv: time=" << simulation_time 
              << ", rirf_steps=" << rirf_steps << std::endl;
#endif

    // Build SystemState from all Chrono bodies (legacy implementation)
    const auto system_state = BuildSystemStateFromChronoBodies(bodies_);

    // Prepare body velocities for recording
    std::vector<std::vector<double>> body_velocities(num_bodies_);
    for (int b = 0; b < num_bodies_; ++b) {
        const auto& body_state = system_state.bodies[b];
        body_velocities[b] = {
            body_state.linear_velocity[0], body_state.linear_velocity[1], body_state.linear_velocity[2],
            body_state.angular_velocity[0], body_state.angular_velocity[1], body_state.angular_velocity[2]
        };
    }

    // Record velocities in the convolution module
    radiation_convolution_->RecordVelocities(simulation_time, body_velocities);

    // Also update legacy history buffers for compatibility (used by GetRIRFval and other legacy code)
    time_history_.insert(time_history_.begin(), simulation_time);
    for (int b = 0; b < num_bodies_; ++b) {
        velocity_history_[b].insert(velocity_history_[b].begin(), body_velocities[b]);
    }

    // Create lambda to get RIRF values (for Baseline mode)
    auto get_rirf_val = [this](int row, int col, int step) -> double {
        return GetRIRFval(row, col, step);
    };

    // Compute forces using the convolution module
    const std::vector<Eigen::Tensor<double, 3>>* processed_rirf = 
        (convolution_mode_ == RadiationConvolutionMode::TaperedDirect && rirf_processed_ready_) 
        ? &rirf_processed_ : nullptr;

    force_radiation_damping_ = radiation_convolution_->ComputeForces(get_rirf_val, processed_rirf);

#ifdef HYDROCHRONO_DEBUG
    std::cout << "[HYDRO_DEBUG] Radiation conv complete: max_force=" 
              << *std::max_element(force_radiation_damping_.begin(), force_radiation_damping_.end()) << std::endl;
#endif

    profile_stats_.radiation_seconds += std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - __t0).count();
    profile_stats_.radiation_calls++;
    return force_radiation_damping_;
}

double TestHydro::GetRIRFval(int row, int col, int st) {
    if (row < 0 || row >= kDofPerBody * num_bodies_ || col < 0 || col >= kDofPerBody * num_bodies_ || st < 0 ||
        st >= file_info_.GetRIRFDims(2)) {
        throw std::out_of_range("rirfval index out of range in TestHydro");
    }

    int body_index = row / kDofPerBody;
    int col_dof    = col % kDofPerBody;
    int row_dof    = row % kDofPerBody;

    if (convolution_mode_ == RadiationConvolutionMode::TaperedDirect) {
        EnsureProcessedRIRF();
        // processed tensor is scaled by rho already
        const auto& tensor = rirf_processed_[body_index];
        return tensor(row_dof, col, st);
    }

    return file_info_.GetRIRFVal(body_index, row_dof, col, st);
}

// ------------------------------------------------------------
// SECTION: Wave excitation
// ------------------------------------------------------------
// Delegates to WaveBase hierarchy for excitation force computation.
// TODO: Extract to ExcitationModel class in future refactor.

Eigen::VectorXd TestHydro::ComputeForceWaves() {
    auto __t0 = std::chrono::steady_clock::now();
    // Ensure bodies_ is not empty
    if (bodies_.empty()) {
        throw std::runtime_error("bodies_ array is empty in ComputeForceWaves");
    }

    force_waves_ = user_waves_->GetForceAtTime(bodies_[0]->GetChTime());

#ifdef HYDROCHRONO_DEBUG
    std::cout << "[HYDRO_DEBUG] Wave excitation: max_force=" 
              << force_waves_.maxCoeff() << std::endl;
#endif

    profile_stats_.waves_seconds += std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - __t0).count();
    profile_stats_.waves_calls++;
    return force_waves_;
}

// ------------------------------------------------------------
// SECTION: Main force evaluation (Chrono callback entry point)
// ------------------------------------------------------------
// Called by Chrono via ForceFunc6d/ComponentFunc callbacks.
// Coordinates all force components and caches results per time step.
// TODO: This will become a thin adapter delegating to HydroSystem.

double TestHydro::CoordinateFuncForBody(int b, int dof_index) {
    if (dof_index < 0 || dof_index >= kDofPerBody || b < 1 || b > num_bodies_) {
        throw std::out_of_range("Invalid index in CoordinateFuncForBody");
    }

    // Adjusting for 1-indexed body number
    const int body_num_offset = kDofPerBody * (b - 1);
    const int total_dofs      = kDofPerBody * num_bodies_;

    // Ensure the bodies_ vector isn't empty and the first element isn't null
    if (bodies_.empty() || !bodies_[0]) {
        throw std::runtime_error("bodies_ array is empty or invalid in CoordinateFuncForBody");
    }

    // Time-step caching: reuse forces if already computed this step
    if (bodies_[0]->GetChTime() == prev_time) {
        return total_force_[body_num_offset + dof_index];
    }

    // Update time and reset forces for this time step
    prev_time = bodies_[0]->GetChTime();
    std::fill(total_force_.begin(), total_force_.end(), 0.0);
    std::fill(force_hydrostatic_.begin(), force_hydrostatic_.end(), 0.0);
    std::fill(force_radiation_damping_.begin(), force_radiation_damping_.end(), 0.0);
    std::fill(force_waves_.begin(), force_waves_.end(), 0.0);

    // Compute all force components (multibody: all bodies computed together)
    force_hydrostatic_       = ComputeForceHydrostatics();
    force_radiation_damping_ = ComputeForceRadiationDampingConv();
    force_waves_             = ComputeForceWaves();

    // Accumulate total force (multibody: sum over all DOFs for all bodies)
    // Note: radiation damping is subtracted (damping opposes motion)
    for (int index = 0; index < total_dofs; index++) {
        total_force_[index] = force_hydrostatic_[index] - force_radiation_damping_[index] + force_waves_[index];
    }

#ifdef HYDROCHRONO_DEBUG
    std::cout << "[HYDRO_DEBUG] CoordinateFuncForBody: body=" << b << ", dof=" << dof_index 
              << ", time=" << prev_time << ", force=" << total_force_[body_num_offset + dof_index] << std::endl;
#endif

    if (body_num_offset + dof_index < 0 || body_num_offset >= total_dofs) {
        throw std::out_of_range("Accessing out-of-bounds index in CoordinateFuncForBody");
    }

    return total_force_[body_num_offset + dof_index];
}