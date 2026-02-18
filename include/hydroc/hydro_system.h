/*********************************************************************
 * @file  hydro_system.h
 * @brief HydroSystem façade: attaches hydrodynamic forces to Chrono bodies.
 *
 * MAIN TYPES:
 *   - HydroSystem: Primary hydrodynamics façade (adapter over HydroForces)
 *   - ForceFunc6d: Wraps 6-DOF force/torque callbacks per body
 *   - ComponentFunc: Per-DOF force function for Chrono's ChForce system
 *
 * ROLE: This is the main entry point for C++ usage. Internally delegates
 * to HydroForces (Chrono-free core) and ChronoHydroCoupler. Used by both
 * C++ tests/demos and the YAML runner (hydrochrono.exe).
 *
 * KEY ASSUMPTIONS: 6 DOF per body, bodies named "body1", "body2", etc.
 *********************************************************************/

#ifndef HYDRO_SYSTEM_H
#define HYDRO_SYSTEM_H

// ─────────────────────────────────────────────────────────────────────────────
// Added Mass Implementation Toggle
// ─────────────────────────────────────────────────────────────────────────────
// By default, HydroChrono uses Chrono's built-in ChLoadHydrodynamics for
// infinite-frequency added mass. To switch to HydroChrono's legacy
// ChLoadAddedMass (ChLoadCustomMultiple-based) implementation, uncomment:
//
// #define HYDROCHRONO_USE_LEGACY_ADDED_MASS
//

// Include hydro_types.h FIRST to ensure BodyForces and GeneralizedForce are available
// before any other includes that might conflict (e.g., config/hydro_config.h)
#include <hydroc/core/hydro_types.h>

// Standard includes
#include <cstdio>
#include <filesystem>
#include <limits>
#include <memory>
#include <vector>

// Eigen includes (for GeneralizedForce and BodyForces)
#include <Eigen/Dense>

// Chrono library includes
// Note: We include only the headers strictly necessary for the types used in this header.
// Implementation-only Chrono headers are included in hydro_system.cpp.
// ChForce.h includes ChFunction internally, so we don't need a separate include.
#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChForce.h>
#ifdef HYDROCHRONO_USE_LEGACY_ADDED_MASS
#include <chrono/physics/ChLoadContainer.h>
#endif

// Hydroc library includes
#include <hydroc/io/h5_reader.h>
#include <hydroc/waves/wave_base.h>
#include <hydroc/waves/regular_wave.h>
#include <hydroc/waves/irregular_wave.h>

// Public hydroc includes (for SystemState)
#include <hydroc/core/system_state.h>

// Radiation types (canonical definitions - single source of truth)
#include <hydroc/radiation/radiation_types.h>

// Force component interface (for radiation component)
#include <hydroc/core/force_component.h>

namespace hydrochrono::hydro {
// Forward declarations
class HydrostaticsComponent;
class RadiationComponent;
class RadiationStateSpaceComponent;
class ExcitationComponent;
class HydroForces;
class ChronoHydroCoupler;
class IHydroForceComponent;
// Types GeneralizedForce and BodyForces are defined in hydro_types.h (included above)
}

class ForceFunc6d;
class HydroSystem;

// ─────────────────────────────────────────────────────────────────────────────
// Degrees of Freedom Constant
// ─────────────────────────────────────────────────────────────────────────────
/// Number of degrees of freedom per rigid body: surge, sway, heave, roll, pitch, yaw.
constexpr int kDofPerBody = 6;

/**
 * @brief Per-DOF force function for Chrono's ChForce system.
 *
 * Inherits from chrono::ChFunction to provide time-dependent force values.
 * Each instance represents a single degree of freedom (surge, sway, heave, roll, pitch, yaw).
 */
class ComponentFunc : public chrono::ChFunction {
  public:
    /**
     * @brief Default constructor.
     * @note Do not use this constructor.
     */
    ComponentFunc();

    /**
     * @brief Copy constructor.
     */
    ComponentFunc(const ComponentFunc& old);

    /**
     * @brief Construct a new ComponentFunc object.
     *
     * @param b Pointer to ForceFunc6d member object.
     * @param i Index for which component this ComponentFunc object refers to. Corresponds
     * to the force's degree of freedom, where:
     * (0,1,2,3,4,5) -> (surge, sway, heave, roll, pitch, yaw).
     */
    ComponentFunc(ForceFunc6d* b, int i);

    /**
     * @brief Clone the object.
     * @note This is a required override since ComponentFunc inherits from ChFunction.
     * @return A cloned ComponentFunc object.
     */
    virtual ComponentFunc* Clone() const override;

    /**
     * @brief Calculate the force value for a specific time.
     *
     * @param x Time from simulation.
     * @return Force on body in the specified degree of freedom at time x.
     */
    virtual double GetVal(double x) const override;

  private:
    ForceFunc6d* base_ = nullptr;  ///< Pointer to the full 6D force on the body.
    int index_ = 0;          ///< Index representing force degree of freedom on the body.
};

/**
 * @brief Organizes the functional (time-dependent) forces in each degree of freedom (6 total) for a body.
 */
class ForceFunc6d {
  public:
    /**
     * @brief Initializes an array of ComponentFunc objects and pointers to each force/torque.
     */
    ForceFunc6d();

    /**
     * @brief Initializes hydro force info from H5FileInfo and the ChBody this force will be applied to.
     *
     * @param object The body in the system to which this 6-dimensional force is being applied.
     * @param all_hydro_forces_user The HydroSystem class where the total force on all bodies is calculated.
     */
    ForceFunc6d(std::shared_ptr<chrono::ChBody> object, HydroSystem* all_hydro_forces_user);

    /**
     * @brief Copy constructor that ensures the force is only added to a body once.
     *
     * @note Avoid using the default copy constructor. ApplyForceAndTorqueToBody() should only ever be applied once.
     *
     * @param old ForceFunc6d object to copy from.
     */
    ForceFunc6d(const ForceFunc6d& old);

    /**
     * @brief Calculates the force on a given degree of freedom.
     *
     * @param i Index corresponding to the degree of freedom. Assumed to be 0-5 only.
     * @return Value of force in the i-th degree of freedom for the associated body.
     */
    double CoordinateFunc(int i);

  private:
    /**
     * @brief Initializes force components.
     */
    void SetForce();

    /**
     * @brief Initializes torque components.
     */
    void SetTorque();

    /**
     * @brief Adds this force to the body's list of applied forces.
     *
     * @warning Ensure this function is not called multiple times for the same force.
     */
    void ApplyForceAndTorqueToBody();

    std::shared_ptr<chrono::ChBody> body_;          ///< Pointer to the body this force is applied to.
    
    /// Body index in the system (1-indexed, parsed from Chrono body name "bodyN").
    /// This matches Chrono's body naming convention where bodies are named body1, body2, etc.
    int b_num_ = 0;
    
    ComponentFunc forces_[kDofPerBody];                       ///< Forces for each degree of freedom.
    std::shared_ptr<ComponentFunc> force_ptrs_[kDofPerBody];  ///< Pointers to the forces.
    std::shared_ptr<chrono::ChForce> chrono_force_;  ///< Chrono force for the body.
    std::shared_ptr<chrono::ChForce> chrono_torque_; ///< Chrono torque for the body.
    HydroSystem* all_hydro_forces_ = nullptr;                 ///< Pointer to HydroSystem for calculations.
};

#ifdef HYDROCHRONO_USE_LEGACY_ADDED_MASS
// Forward declaration of ChLoadAddedMass (defined in added_mass.h, included in .cpp)
class ChLoadAddedMass;
#endif

// Lightweight hydrodynamics profiling stats
struct HydroProfileStats {
    double hydrostatics_seconds = 0.0;
    double radiation_seconds    = 0.0;
    double waves_seconds        = 0.0;
    int hydrostatics_calls      = 0;
    int radiation_calls         = 0;
    int waves_calls             = 0;
};

// ═══════════════════════════════════════════════════════════════════════════════
//
//  HydroSystem — Primary Hydrodynamics Façade / Adapter
//
// ═══════════════════════════════════════════════════════════════════════════════
/**
 * @brief Primary hydrodynamics façade for HydroChrono applications.
 *
 * HydroSystem is the main user-facing object for attaching hydrodynamic forces to
 * Chrono bodies. This is the production hydrodynamics adapter used by all workflows.
 *
 * @note The historical names "TestHydro" and "HydroForces" are retained as aliases
 *       for backward compatibility. New code should use `HydroSystem`.
 *
 * ## What It Does
 *
 * HydroSystem bridges three layers:
 *   1. **Chrono bodies** — The physical bodies in the Chrono simulation.
 *   2. **HydroForces (core)** — Chrono-free hydrodynamic force computation.
 *   3. **Chrono force callbacks** — ChForce functions that Chrono calls each timestep.
 *
 * ## Typical Usage
 *
 * @code{.cpp}
 *   // Create Chrono bodies
 *   auto body1 = chrono_types::make_shared<ChBody>();
 *   body1->SetName("body1");
 *   system.Add(body1);
 *
 *   // Attach hydrodynamic forces (with or without waves)
 *   std::vector<std::shared_ptr<ChBody>> bodies = {body1};
 *   HydroSystem hydro(bodies, "path/to/hydro_data.h5");
 *   hydro.AddWaves(std::make_shared<RegularWave>(...));
 *
 *   // Run simulation — Chrono automatically queries HydroSystem for forces
 * @endcode
 *
 * ## Used By
 *
 *   - **C++ tests and demos** — Directly instantiate HydroSystem with bodies and H5 file.
 *   - **YAML runner (hydrochrono.exe)** — SetupHydroFromYAML creates HydroSystem internally.
 *
 * ## Architecture Notes
 *
 *   - Internally delegates to HydroForces + ChronoHydroCoupler for force computation.
 *   - Maintains velocity history for radiation damping convolution.
 *   - Caches forces per timestep to avoid redundant computation.
 *
 * ## Body Indexing Convention
 *
 *   Chrono bodies must be named "body1", "body2", etc. The numeric suffix is parsed
 *   to determine body order. Internally, body indices are **1-based** in the callback
 *   interface (CoordinateFuncForBody) to match this naming convention.
 *
 * @see TestHydro (backwards-compatibility alias)
 * @see HydroForces (Chrono-free hydrodynamic core in hydrochrono::hydro namespace)
 * @see ChronoHydroCoupler (extracts state from Chrono bodies)
 */
class HydroSystem {
  public:
    HydroSystem() = delete;

    /**
     * @brief Main constructor for initializing the HydroSystem class.
     *
     * Sets up vector of bodies, h5 file info, and hydro inputs. If no waves are given,
     * this constructor defaults to using NoWave.
     *
     * @param user_bodies List of pointers to bodies for the hydro forces.
     * @param h5_file_name Name of the h5 file where hydro data is stored.
     * @param waves WaveBase object. Defaults to NoWave if not provided.
     */
    HydroSystem(std::vector<std::shared_ptr<chrono::ChBody>> user_bodies,
                std::string h5_file_name,
                std::shared_ptr<WaveBase> waves = std::make_shared<NoWave>());

    // Destructor (defined in .cpp to allow unique_ptr to incomplete type)
    ~HydroSystem();

    // Deleted copy constructor and assignment operator for safety.
    HydroSystem(const HydroSystem& old) = delete;
    HydroSystem& operator=(const HydroSystem& rhs) = delete;

    /**
     * @brief Adds waves class to force calculations depending on if regular or irregular waves.
     *
     * Also initializes h5 data for wave force class.
     *
     * @param waves the specific WaveBase class to add to the system
     */
    void AddWaves(std::shared_ptr<WaveBase> waves);

    /**
     * @brief Legacy API: Computes hydrostatic forces for a 6N dimensional system.
     *
     * NOTE: Retained for backward compatibility. The main force path now uses HydroForces
     * (via CoordinateFuncForBody). This method is not called during normal simulation.
     *
     * @return 6N dimensional force for 6 DOF and N bodies in system.
     */
    std::vector<double> ComputeForceHydrostatics();

    /**
     * @brief Legacy API: Computes radiation damping forces for a 6N dimensional system.
     *
     * NOTE: Retained for backward compatibility. The main force path now uses HydroForces
     * (via CoordinateFuncForBody). This method is not called during normal simulation.
     *
     * Returns positive damping magnitudes (legacy convention); the main path applies
     * the correct sign internally via RadiationComponent.
     *
     * @return 6N dimensional force for 6 DOF and N bodies in system.
     */
    std::vector<double> ComputeForceRadiationDampingConv();

    /**
     * @brief Legacy API: Computes wave excitation forces for a 6N dimensional system.
     *
     * NOTE: Retained for backward compatibility. The main force path now uses HydroForces
     * (via CoordinateFuncForBody). This method is not called during normal simulation.
     *
     * @return 6N dimensional force for 6 DOF and N bodies in system.
     */
    Eigen::VectorXd ComputeForceWaves();
    // Expose the wave object (read-only) so callers can query inputs if needed
    std::shared_ptr<WaveBase> GetWave() const { return user_waves_; }

    /**
     * @brief Fetches the RIRF value from the h5 file based on the provided indices.
     *
     * @param row Index representing the body number and DOF index [0,...,5,...6N-1] for rows of RIRF.
     * @param col Column index in RIRF matrix [0,...,5,...6N-1].
     * @param st Index representing the timestep in RIRF, usually in the range [0,...1000].
     *
     * @return The value from the h5 file RIRF matrix.
     */
    double GetRIRFval(int row, int col, int st);

    // ─────────────────────────────────────────────────────────────────────────
    // Radiation Configuration
    // ─────────────────────────────────────────────────────────────────────────
    // These types are defined in the radiation module (single source of truth):
    //   - hydrochrono::hydro::RadiationMethod (enum: kRirfConvolution, kStateSpace)
    //   - hydrochrono::hydro::RadiationConvolutionMode (enum: Baseline, TaperedDirect)
    //   - hydrochrono::hydro::TaperedDirectOptions (struct)
    //   - hydrochrono::hydro::StateSpaceOptions (struct)
    // See: include/hydroc/radiation/radiation_types.h

    /**
     * @brief Set the radiation damping calculation method.
     * 
     * Selects the overall approach:
     *   - kRirfConvolution: Direct RIRF convolution (default)
     *   - kStateSpace: State-space exponential approximation
     * 
     * @note State-space method not yet implemented. This value is stored but
     *       currently has no effect; runtime always uses RIRF convolution.
     * 
     * @param method RadiationMethod::kRirfConvolution or ::kStateSpace
     */
    void SetRadiationMethod(hydrochrono::hydro::RadiationMethod method);

    /**
     * @brief Set options for state-space radiation approximation.
     * 
     * Only applies when RadiationMethod == kStateSpace.
     * 
     * @note State-space method not yet implemented. Options stored for future use.
     * 
     * @param opts StateSpaceOptions struct with max_order and r2_threshold
     */
    void SetStateSpaceOptions(const hydrochrono::hydro::StateSpaceOptions& opts);

    /**
     * @brief Set the radiation convolution mode. Default is Baseline.
     * 
     * Only applies when RadiationMethod == kRirfConvolution.
     * 
     * @param mode RadiationConvolutionMode::Baseline or ::TaperedDirect
     */
    void SetRadiationConvolutionMode(hydrochrono::hydro::RadiationConvolutionMode mode);

    /**
     * @brief Set options for TaperedDirect preprocessing.
     * 
     * Only applies when RadiationMethod == kRirfConvolution and
     * RadiationConvolutionMode == TaperedDirect.
     * 
     * @param opts TaperedDirectOptions struct with smoothing, tapering, and export settings
     */
    void SetTaperedDirectOptions(const hydrochrono::hydro::TaperedDirectOptions& opts);

    /**
     * @brief Set the directory where diagnostics (e.g., CSVs) should be written.
     */
    void SetDiagnosticsOutputDirectory(const std::string& dir) {
        diagnostics_output_dir_ = dir;
        InvalidateRadiationComponent();  // Invalidate component to recreate with new settings
    }

    /**
     * @brief Calculates or retrieves the total force on a specific body in a particular degree of freedom.
     *
     * This is the main entry point called by Chrono's ChForce callbacks during simulation.
     * Forces are computed once per timestep and cached; subsequent calls within the same
     * timestep return the cached value.
     *
     * @note **1-Based Body Indexing**: The body index `b` is 1-based (body1 → b=1, body2 → b=2, etc.)
     *       to match Chrono's body naming convention. This is intentional and matches how
     *       ForceFunc6d parses body names like "body1", "body2", etc.
     *
     * @param b Body index (1-based: valid range is [1, num_bodies]).
     * @param i Degree of Freedom (DOF) index (0-based: 0=surge, 1=sway, 2=heave, 3=roll, 4=pitch, 5=yaw).
     *
     * @return Force component for body 'b' in DOF 'i' (Newtons for forces, Newton-meters for torques).
     *
     * @throws std::out_of_range if b or i are out of valid range.
     */
    double CoordinateFuncForBody(int b, int i);

    // Hydrodynamics profiling accessors
    HydroProfileStats GetProfileStats() const { return profile_stats_; }

    /**
     * @brief Enable or disable profiling in HydroForces.
     * 
     * When enabled, HydroForces measures timing for each force component.
     * Disabled by default to avoid overhead in normal runs.
     * 
     * @param enabled True to enable profiling, false to disable
     */
    void SetProfilingEnabled(bool enabled);

    // Compare mode: legacy debugging feature, retained for API compatibility.
    // No longer affects behavior since the main path now uses HydroForces exclusively.
    void SetCompareMode(bool enable) { compare_mode_ = enable; }
    bool GetCompareMode() const { return compare_mode_; }

  private:
    // Class properties related to the body and hydrodynamics
    std::vector<std::shared_ptr<chrono::ChBody>> bodies_;
    int num_bodies_;
    HydroData file_info_;
    std::vector<ForceFunc6d> force_per_body_;
    std::shared_ptr<WaveBase> user_waves_;

    // Force components vectors
    std::vector<double> force_hydrostatic_;
    std::vector<double> force_radiation_damping_;
    Eigen::VectorXd force_waves_;
    std::vector<double> total_force_;  // Saved force per timestep to reduce redundant calculations

    // Additional properties related to equilibrium and hydrodynamics
    std::vector<double> equilibrium_;
    std::vector<double> cb_minus_cg_;
    Eigen::VectorXd rirf_time_vector;  // Assumed consistent for each body
    Eigen::VectorXd rirf_width_vector;

    // Time tracking for force caching
    double prev_time;
    
    // Debug counter for printing first few timesteps
    int debug_call_count_ = 0;
    
    // Track if we've already reported NaN (to avoid spamming)
    bool nan_reported_ = false;

    // Cached SystemState: built once per time step and reused by all force computations
    hydrochrono::hydro::SystemState cached_state_;
    double cached_state_time_ = std::numeric_limits<double>::quiet_NaN();

#ifdef HYDROCHRONO_USE_LEGACY_ADDED_MASS
    // Added mass related properties (legacy ChLoadAddedMass path)
    std::shared_ptr<chrono::ChLoadContainer> my_loadcontainer;
    std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
#endif

    /**
     * @brief Returns the cached SystemState for the given time.
     *
     * If cached_state_time_ matches the requested time, returns the cached state.
     * Otherwise, builds a fresh state and updates the cache (should not happen
     * in normal flow since CoordinateFuncForBody builds it first).
     *
     * @param time Simulation time to get state for
     * @return Reference to the cached SystemState
     */
    const hydrochrono::hydro::SystemState& GetCachedSystemState(double time);

    // Hydrodynamics profiling data (accumulated over run)
    HydroProfileStats profile_stats_;

    // Hydrostatics force component
    std::unique_ptr<hydrochrono::hydro::HydrostaticsComponent> hydrostatics_component_;

    // Radiation damping force component (convolution or state-space based on radiation_method_)
    std::unique_ptr<hydrochrono::hydro::IHydroForceComponent> radiation_component_;

    // Wave excitation force component
    std::unique_ptr<hydrochrono::hydro::ExcitationComponent> excitation_component_;

    // HydroForces + ChronoHydroCoupler (persistent, constructed once)
    std::unique_ptr<hydrochrono::hydro::HydroForces> hydro_forces_;
    std::unique_ptr<hydrochrono::hydro::ChronoHydroCoupler> chrono_coupler_;

    // Legacy compare mode flag: retained for API compatibility; no longer used.
    bool compare_mode_ = false;

    // Profiling enable flag: stored here because chrono_coupler_ is created lazily.
    bool profiling_enabled_ = false;

    // ─────────────────────────────────────────────────────────────────────────
    // Radiation configuration (uses canonical types from radiation module)
    // ─────────────────────────────────────────────────────────────────────────
    
    // Top-level method selection (state-space not yet implemented)
    hydrochrono::hydro::RadiationMethod radiation_method_ = 
        hydrochrono::hydro::RadiationMethod::kRirfConvolution;
    hydrochrono::hydro::StateSpaceOptions state_space_opts_;
    
    // Convolution kernel preprocessing (only applies when method == kRirfConvolution)
    hydrochrono::hydro::RadiationConvolutionMode convolution_mode_;
    bool rirf_processed_ready_ = false;
    std::vector<Eigen::Tensor<double, 3>> rirf_processed_; // per body [dof x col x step]
    hydrochrono::hydro::TaperedDirectOptions tapered_opts_;
    std::string diagnostics_output_dir_;

    void EnsureProcessedRIRF();

    // Helper to ensure radiation component exists with current settings
    void EnsureRadiationComponent();
    // Helper to invalidate radiation component (requires full type definition)
    void InvalidateRadiationComponent();
    // Helper to ensure excitation component exists
    void EnsureExcitationComponent();

    // Factory: creates ExcitationComponent with current wave configuration.
    // Used by both EnsureExcitationComponent() and EnsureHydroForcesAndCoupler()
    // to ensure consistent construction.
    std::unique_ptr<hydrochrono::hydro::ExcitationComponent> CreateExcitationComponent() const;

    // Factory: creates HydrostaticsComponent with current equilibrium/stiffness data.
    // Used by both the HydroSystem constructor and EnsureHydroForcesAndCoupler()
    // to ensure consistent construction.
    std::unique_ptr<hydrochrono::hydro::HydrostaticsComponent> CreateHydrostaticsComponent() const;

    // Factory: creates radiation force component based on radiation_method_.
    // Returns RadiationComponent (convolution) or RadiationStateSpaceComponent (state-space).
    // Used by both EnsureRadiationComponent() and EnsureHydroForcesAndCoupler()
    // to ensure consistent construction.
    // Note: Each instance owns its own internal state (they are independent).
    std::unique_ptr<hydrochrono::hydro::IHydroForceComponent> CreateRadiationComponent() const;

    // Internal helper: constructs hydro_forces_ and chrono_coupler_ once.
    // Subsequent calls are no-ops. Called automatically by CoordinateFuncForBody().
    void EnsureHydroForcesAndCoupler();
};

// ─────────────────────────────────────────────────────────────────────────────
// Backwards-Compatibility Aliases
// ─────────────────────────────────────────────────────────────────────────────
/**
 * @brief Backwards-compatibility alias for HydroSystem.
 *
 * TestHydro is the historical name for the main hydrodynamics façade.
 * It is retained as an alias for backward compatibility with existing code.
 *
 * @note Prefer `HydroSystem` for new code.
 *
 * @deprecated Use HydroSystem instead.
 */
using TestHydro = HydroSystem;

/**
 * @brief Backwards-compatibility alias for HydroSystem.
 *
 * HydroForces was the previous name for the main hydrodynamics façade
 * (before the HydroSystem/HydroForces rename). It is retained as an alias
 * for backward compatibility with existing code.
 *
 * @note Prefer `HydroSystem` for new code.
 * @note Not to be confused with hydrochrono::hydro::HydroForces, which is the
 *       internal Chrono-free force computation engine.
 *
 * @deprecated Use HydroSystem instead.
 */
using HydroForces = HydroSystem;

#endif

