/*********************************************************************
 * @file  hydro_forces.cpp
 *
 * @brief Implementation of TestHydro main class and helper classes
 * ComponentFunc and ForceFunc6d.
 *
 * ARCHITECTURE:
 * TestHydro is a thin adapter over HydroSystem + ChronoHydroCoupler.
 *
 * Force computation flow:
 *   1. Chrono calls CoordinateFuncForBody(body, dof) via ChForce callbacks.
 *   2. CoordinateFuncForBody detects new timesteps and calls chrono_coupler_->Evaluate().
 *   3. ChronoHydroCoupler extracts system state from Chrono bodies and invokes HydroSystem.
 *   4. HydroSystem evaluates all force components and returns BodyForces.
 *   5. Forces are flattened into the legacy total_force_ array for Chrono consumption.
 *
 * Sign convention: total = hydrostatics - radiation + waves
 *   (RadiationComponent applies the negative sign internally since damping opposes motion.)
 *
 * COMPONENT CONSTRUCTION:
 * Force components are built via factory methods (single source of truth):
 *   - CreateHydrostaticsComponent()
 *   - CreateRadiationComponent()
 *   - CreateExcitationComponent()
 *
 * LEGACY API:
 * The following methods are retained for backward compatibility but are NOT used
 * by the main force path (CoordinateFuncForBody):
 *   - ComputeForceHydrostatics()
 *   - ComputeForceRadiationDampingConv()
 *   - ComputeForceWaves()
 *
 * MAIN RESPONSIBILITIES:
 * - TestHydro: Façade over HydroSystem; provides Chrono force callbacks
 * - ForceFunc6d: Wraps Chrono ChForce/ChTorque callbacks for each body
 * - ComponentFunc: Provides per-DOF force values to Chrono's force system
 *
 * INTERACTIONS:
 * - Reads HDF5 hydrodynamic data via H5FileInfo (equilibrium, stiffness, RIRF kernels)
 * - Uses WaveBase hierarchy for wave excitation forces
 * - Applies forces to Chrono bodies through ChForce/ChLoadAddedMass
 * - RadiationComponent maintains velocity history for radiation convolution
 *
 * KEY ASSUMPTIONS:
 * - All bodies share the same H5 file (multibody data in single file)
 * - Bodies are 1-indexed in ForceFunc6d interface (legacy)
 * - 6 DOF per body (surge, sway, heave, roll, pitch, yaw)
 * - Forces computed once per time step, cached via prev_time check
 *
 * DEBUG INSTRUMENTATION:
 * To enable debug output, compile with -DHYDROCHRONO_DEBUG.
 *********************************************************************/

#include "hydroc/hydro_forces.h"
#include <hydroc/coupling/added_mass.h>
#include <hydroc/io/h5_reader.h>
#include <hydroc/waves/wave_base.h>
#include <hydroc/waves/regular_wave.h>
#include <hydroc/waves/irregular_wave.h>
#include <hydroc/logging.h>
#include <hydroc/core/system_state.h>
#include "hydro/core/chrono_state_utils.h"
#include "hydro/force_components/hydrostatics_component.h"
#include "hydro/force_components/radiation_component.h"
#include "hydro/force_components/excitation_component.h"
#include <hydroc/core/hydro_system.h>
#include <hydroc/coupling/chrono_coupler.h>
#include "hydro/radiation/radiation_rirf_processing.h"

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
// TODO: This class will be refactored into modular force components.
// Current structure mixes all physics components together.

// Explicit destructor definition (needed for unique_ptr to incomplete type)
TestHydro::~TestHydro() = default;

TestHydro::TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
                     std::string h5_file_name,
                     std::shared_ptr<WaveBase> waves)
    : bodies_(user_bodies),
      num_bodies_(bodies_.size()),
      file_info_(H5FileInfo(h5_file_name, num_bodies_).ReadH5Data()),
      hydro_system_(nullptr),
      chrono_coupler_(nullptr) {
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

    // Initialize force vectors
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

    // Initialize hydrostatics component (uses shared factory for consistent construction)
    hydrostatics_component_ = CreateHydrostaticsComponent();

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
// SECTION: Cached SystemState helper
// ------------------------------------------------------------
// Returns the cached SystemState for the given time. In normal flow,
// CoordinateFuncForBody builds the state at the start of each timestep.
// This helper provides a fallback if called out of order.

const hydrochrono::hydro::SystemState& TestHydro::GetCachedSystemState(double time) {
    // Check if cached state is valid for this time
    if (cached_state_time_ == time) {
        return cached_state_;
    }

    // Fallback: build state now (should not happen in normal flow)
    hydrochrono::hydro::BuildSystemStateFromChronoBodies(bodies_, cached_state_);
    cached_state_time_ = time;
    return cached_state_;
}

// ------------------------------------------------------------
// SECTION: Hydrostatics (legacy path, delegates to HydrostaticsComponent)
// ------------------------------------------------------------
// NOTE: This method is kept for backward compatibility but is no longer called by
// CoordinateFuncForBody(). The main force path now goes through HydroSystem.

std::vector<double> TestHydro::ComputeForceHydrostatics() {
    auto __t0 = std::chrono::steady_clock::now();
    assert(num_bodies_ > 0);

#ifdef HYDROCHRONO_DEBUG
    const double rho = file_info_.GetRhoVal();
    std::cout << "[HYDRO_DEBUG] ComputeForceHydrostatics: num_bodies=" << num_bodies_ << ", rho=" << rho << std::endl;
#endif

    // Get current simulation time and cached state
    const double simulation_time = bodies_[0]->GetChTime();
    const hydrochrono::hydro::SystemState& system_state = GetCachedSystemState(simulation_time);

    // Compute forces using the hydrostatics component
    hydrochrono::hydro::BodyForces body_forces(num_bodies_);
    for (int b = 0; b < num_bodies_; ++b) {
        body_forces[b].resize(kDofPerBody);
        body_forces[b].setZero();
    }
    hydrostatics_component_->Compute(system_state, simulation_time, body_forces);

    // Convert BodyForces back to legacy flat 6N vector format
    force_hydrostatic_.assign(kDofPerBody * num_bodies_, 0.0);
    for (int b = 0; b < num_bodies_; ++b) {
        const int body_offset = kDofPerBody * b;
        for (int i = 0; i < kDofPerBody; ++i) {
            force_hydrostatic_[body_offset + i] = body_forces[b][i];
        }
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
// SECTION: Radiation convolution (delegates to RadiationComponent)
// ------------------------------------------------------------
// Delegates to RadiationComponent for force computation.
// Legacy history buffers maintained for compatibility.

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

void TestHydro::InvalidateRadiationComponent() {
    radiation_component_.reset();
}

void TestHydro::EnsureRadiationComponent() {
    if (radiation_component_) {
        return;  // Already created
    }

    radiation_component_ = CreateRadiationComponent();
}

void TestHydro::EnsureExcitationComponent() {
    if (excitation_component_) {
        return;  // Already created
    }

    excitation_component_ = CreateExcitationComponent();
}

std::unique_ptr<hydrochrono::hydro::ExcitationComponent> TestHydro::CreateExcitationComponent() const {
    // Single source of truth for ExcitationComponent construction.
    // Both EnsureExcitationComponent() and EnsureHydroSystemAndCoupler() use this.
    return std::make_unique<hydrochrono::hydro::ExcitationComponent>(user_waves_, num_bodies_);
}

std::unique_ptr<hydrochrono::hydro::HydrostaticsComponent> TestHydro::CreateHydrostaticsComponent() const {
    // Single source of truth for HydrostaticsComponent construction.
    // Used by TestHydro constructor and EnsureHydroSystemAndCoupler().
    const auto gravitational_acceleration_ch = bodies_[0]->GetSystem()->GetGravitationalAcceleration();
    const Eigen::Vector3d gravitational_acceleration(
        gravitational_acceleration_ch.x(),
        gravitational_acceleration_ch.y(),
        gravitational_acceleration_ch.z()
    );
    return std::make_unique<hydrochrono::hydro::HydrostaticsComponent>(
        file_info_, num_bodies_, equilibrium_, cb_minus_cg_, gravitational_acceleration);
}

std::unique_ptr<hydrochrono::hydro::RadiationComponent> TestHydro::CreateRadiationComponent() const {
    // Single source of truth for RadiationComponent construction.
    // Used by EnsureRadiationComponent() and EnsureHydroSystemAndCoupler().
    // Each RadiationComponent instance owns its own velocity history.
    //
    // Configuration inputs:
    //   - file_info_: BEM data including RIRF kernels
    //   - num_bodies_: number of bodies in system
    //   - rirf_time_vector, rirf_width_vector: time discretization for RIRF
    //   - convolution_mode_: Baseline or TaperedDirect
    //   - tapered_opts_: preprocessing options for TaperedDirect mode
    //   - diagnostics_output_dir_: where to write debug CSVs

    const int rirf_steps = file_info_.GetRIRFDims(2);

    // Convert TestHydro::RadiationConvolutionMode to hydrochrono::hydro::RadiationConvolutionMode
    hydrochrono::hydro::RadiationConvolutionMode component_mode;
    if (convolution_mode_ == RadiationConvolutionMode::TaperedDirect) {
        component_mode = hydrochrono::hydro::RadiationConvolutionMode::TaperedDirect;
    } else {
        component_mode = hydrochrono::hydro::RadiationConvolutionMode::Baseline;
    }

    // Convert TestHydro::TaperedDirectOptions to hydrochrono::hydro::TaperedDirectOptions
    hydrochrono::hydro::TaperedDirectOptions component_opts;
    component_opts.smoothing = tapered_opts_.smoothing;
    component_opts.window_length = tapered_opts_.window_length;
    component_opts.rirf_end_time = tapered_opts_.rirf_end_time;
    component_opts.taper_start_percent = tapered_opts_.taper_start_percent;
    component_opts.taper_end_percent = tapered_opts_.taper_end_percent;
    component_opts.taper_final_amplitude = tapered_opts_.taper_final_amplitude;
    component_opts.export_plot_csv = tapered_opts_.export_plot_csv;

    return std::make_unique<hydrochrono::hydro::RadiationComponent>(
        file_info_, num_bodies_, rirf_steps, rirf_time_vector, rirf_width_vector,
        component_mode, component_opts, diagnostics_output_dir_);
}

// ------------------------------------------------------------
// Internal helpers for HydroSystem + ChronoHydroCoupler path
// ------------------------------------------------------------
// HydroSystem + ChronoHydroCoupler are persistent (constructed once).
// The HydroSystem owns its own force components with independent state.
// This path is used by the comparison harness and will eventually replace
// the legacy per-component methods.

void TestHydro::EnsureHydroSystemAndCoupler() {
    if (hydro_system_ && chrono_coupler_) {
        return;  // Already constructed; do not recreate
    }

    // Build force components for HydroSystem (independent from legacy components)
    std::vector<std::unique_ptr<hydrochrono::hydro::IHydroForceComponent>> components;

    // Hydrostatics component (uses shared factory for consistent construction)
    components.push_back(CreateHydrostaticsComponent());

    // Radiation component (uses shared factory for consistent construction)
    components.push_back(CreateRadiationComponent());

    // Excitation component (uses shared factory for consistent construction)
    components.push_back(CreateExcitationComponent());

    // Construct HydroSystem (takes ownership of components)
    hydro_system_ = std::make_unique<hydrochrono::hydro::HydroSystem>(
        num_bodies_, std::move(components));

    // Create shared_ptr alias for ChronoHydroCoupler (empty deleter - unique_ptr owns lifetime)
    std::shared_ptr<hydrochrono::hydro::HydroSystem> hydro_system_shared(
        hydro_system_.get(), [](hydrochrono::hydro::HydroSystem*) {
            // Empty deleter - unique_ptr owns the lifetime
        });

    // Construct ChronoHydroCoupler
    chrono_coupler_ = std::make_unique<hydrochrono::hydro::ChronoHydroCoupler>(
        hydro_system_shared, bodies_);

    // Apply deferred profiling setting (may have been set before coupler was created)
    if (profiling_enabled_) {
        chrono_coupler_->SetProfilingEnabled(true);
    }
}


// Legacy radiation convolution computation.
// RadiationComponent is the single source of truth for radiation history and forces.
// NOTE: This method is kept for backward compatibility but is no longer called by
// CoordinateFuncForBody(). The main force path now goes through HydroSystem.
std::vector<double> TestHydro::ComputeForceRadiationDampingConv() {
    auto __t0 = std::chrono::steady_clock::now();
    const int total_dofs = kDofPerBody * num_bodies_;

    assert(total_dofs > 0);

    // Ensure radiation component exists with current settings
    EnsureRadiationComponent();

    // Current time and cached state
    const double simulation_time = bodies_[0]->GetChTime();
    const hydrochrono::hydro::SystemState& system_state = GetCachedSystemState(simulation_time);

#ifdef HYDROCHRONO_DEBUG
    const int rirf_steps = file_info_.GetRIRFDims(2);
    std::cout << "[HYDRO_DEBUG] Radiation conv: time=" << simulation_time 
              << ", rirf_steps=" << rirf_steps << std::endl;
#endif

    // RadiationComponent::Compute() handles:
    //   - Recording current velocities into its internal history
    //   - Performing RIRF convolution over the velocity history
    //   - Returning the radiation damping forces (negative, since damping opposes motion)
    hydrochrono::hydro::BodyForces body_forces(num_bodies_);
    for (int b = 0; b < num_bodies_; ++b) {
        body_forces[b].resize(kDofPerBody);
        body_forces[b].setZero();
    }
    radiation_component_->Compute(system_state, simulation_time, body_forces);

    // Convert BodyForces to legacy flat 6N vector format.
    // Negate to return positive damping magnitude (legacy convention: CoordinateFuncForBody
    // used to subtract this, but now HydroSystem handles the sign internally).
    force_radiation_damping_.assign(kDofPerBody * num_bodies_, 0.0);
    for (int b = 0; b < num_bodies_; ++b) {
        const int body_offset = kDofPerBody * b;
        for (int i = 0; i < kDofPerBody; ++i) {
            force_radiation_damping_[body_offset + i] = -body_forces[b][i];  // Negate for legacy compatibility
        }
    }

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
// SECTION: Wave excitation (legacy path, delegates to ExcitationComponent)
// ------------------------------------------------------------
// NOTE: This method is kept for backward compatibility but is no longer called by
// CoordinateFuncForBody(). The main force path now goes through HydroSystem.

Eigen::VectorXd TestHydro::ComputeForceWaves() {
    auto __t0 = std::chrono::steady_clock::now();
    // Ensure bodies_ is not empty
    if (bodies_.empty()) {
        throw std::runtime_error("bodies_ array is empty in ComputeForceWaves");
    }

    // Ensure excitation component exists
    EnsureExcitationComponent();

    // Get current simulation time and cached state
    const double simulation_time = bodies_[0]->GetChTime();
    const hydrochrono::hydro::SystemState& system_state = GetCachedSystemState(simulation_time);

    // Compute forces using the excitation component
    hydrochrono::hydro::BodyForces body_forces(num_bodies_);
    for (int b = 0; b < num_bodies_; ++b) {
        body_forces[b].resize(kDofPerBody);
        body_forces[b].setZero();
    }
    excitation_component_->Compute(system_state, simulation_time, body_forces);

    // Convert BodyForces back to legacy flat Eigen::VectorXd format (6N)
    force_waves_.resize(kDofPerBody * num_bodies_);
    for (int b = 0; b < num_bodies_; ++b) {
        const int body_offset = kDofPerBody * b;
        for (int i = 0; i < kDofPerBody; ++i) {
            force_waves_[body_offset + i] = body_forces[b][i];
        }
    }

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
// Routes force computation through HydroSystem + ChronoHydroCoupler.
// Forces are cached per timestep via prev_time check.

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

    // Update time for this new timestep
    prev_time = bodies_[0]->GetChTime();

    // Build SystemState once for this timestep
    hydrochrono::hydro::BuildSystemStateFromChronoBodies(bodies_, cached_state_);
    cached_state_time_ = prev_time;

    // Ensure HydroSystem + ChronoHydroCoupler are initialized
    EnsureHydroSystemAndCoupler();

    // Compute total forces via HydroSystem.
    // HydroSystem::Evaluate() combines all components:
    //   - Hydrostatics: added with + sign
    //   - Radiation: added with - sign (damping opposes motion, handled inside RadiationComponent)
    //   - Excitation: added with + sign
    // Result: total = hydrostatics - radiation + waves
    hydrochrono::hydro::BodyForces body_forces = chrono_coupler_->Evaluate(prev_time);

    // Copy profiling stats from HydroSystem to TestHydro's profile_stats_
    auto hydro_stats = chrono_coupler_->GetProfileStats();
    profile_stats_.hydrostatics_seconds = hydro_stats.hydrostatics_seconds;
    profile_stats_.hydrostatics_calls   = hydro_stats.hydrostatics_calls;
    profile_stats_.radiation_seconds    = hydro_stats.radiation_seconds;
    profile_stats_.radiation_calls      = hydro_stats.radiation_calls;
    profile_stats_.waves_seconds        = hydro_stats.excitation_seconds;  // excitation == waves
    profile_stats_.waves_calls          = hydro_stats.excitation_calls;

    // Flatten BodyForces into total_force_ (legacy flat 6N format)
    std::fill(total_force_.begin(), total_force_.end(), 0.0);
    for (int body_idx = 0; body_idx < num_bodies_; ++body_idx) {
        const int offset = kDofPerBody * body_idx;
        for (int dof = 0; dof < kDofPerBody; ++dof) {
            total_force_[offset + dof] = body_forces[body_idx][dof];
        }
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

void TestHydro::SetProfilingEnabled(bool enabled) {
    profiling_enabled_ = enabled;
    if (chrono_coupler_) {
        chrono_coupler_->SetProfilingEnabled(enabled);
    }
}