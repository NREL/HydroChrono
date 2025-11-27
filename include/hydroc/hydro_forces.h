#ifndef HYDRO_FORCES_H
#define HYDRO_FORCES_H

/*********************************************************************
 * @file  hydro_forces.h
 *
 * @brief Header file of TestHydro main class and helper classes
 * ComponentFunc and ForceFunc6d.
 *
 * ARCHITECTURE:
 * TestHydro is a thin adapter over HydroSystem + ChronoHydroCoupler.
 * Force computation is handled internally by HydroSystem, which owns:
 *   - HydrostaticsComponent (restoring forces + buoyancy)
 *   - RadiationComponent (RIRF convolution, owns velocity history)
 *   - ExcitationComponent (wave forcing)
 *
 * The sign convention is: total = hydrostatics - radiation + waves.
 * (RadiationComponent applies the negative sign internally since damping
 * opposes motion.)
 *
 * MAIN RESPONSIBILITIES:
 * - TestHydro: Façade over HydroSystem; provides Chrono force callbacks
 * - ForceFunc6d: Wraps 6-DOF force/torque callbacks for Chrono bodies
 * - ComponentFunc: Per-DOF force function for Chrono's ChForce system
 *
 * COMPONENT CONSTRUCTION:
 * Force components are constructed via factory methods:
 *   - CreateHydrostaticsComponent()
 *   - CreateRadiationComponent()
 *   - CreateExcitationComponent()
 * These are the single source of truth for component configuration.
 *
 * INTERACTIONS:
 * - Used by setup_hydro_from_yaml to create and configure TestHydro instances
 * - Called by Chrono during simulation via ChForce callbacks
 * - Reads HDF5 data through H5FileInfo (not directly exposed here)
 * - Uses WaveBase hierarchy for wave excitation (passed in constructor)
 *
 * KEY ASSUMPTIONS:
 * - Bodies are 1-indexed in CoordinateFuncForBody (legacy, from ForceFunc6d)
 * - All bodies share same H5 file (multibody data in single file)
 * - 6 DOF per body (surge, sway, heave, roll, pitch, yaw)
 * - Forces computed once per time step and cached via prev_time check
 *********************************************************************/

// Include hydro_types.h FIRST to ensure BodyForces and GeneralizedForce are available
// before any other includes that might conflict (e.g., hydro_config_types.h)
#include <hydroc/hydro_types.h>

// Standard includes
#include <cstdio>
#include <filesystem>
#include <limits>
#include <memory>
#include <vector>

// Eigen includes (for GeneralizedForce and BodyForces)
#include <Eigen/Dense>

// Chrono library includes
#include <chrono/solver/ChIterativeSolverLS.h>
#include <chrono/solver/ChSolverPMINRES.h>
#include <chrono/timestepper/ChTimestepper.h>
#include <chrono/timestepper/ChTimestepperHHT.h>

#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChForce.h>
#include <chrono/physics/ChLoad.h>
#include <chrono/physics/ChLoadContainer.h>
#include <chrono/physics/ChLoadsBody.h>
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/physics/ChSystemSMC.h>

// Chrono FEA includes
#include <chrono/fea/ChMeshFileLoader.h>

// Hydroc library includes
#include <hydroc/h5fileinfo.h>
#include <hydroc/wave_types.h>
// Note: hydro_types.h is included at the top to avoid header guard conflicts

// Public hydroc includes (for SystemState)
#include <hydroc/system_state.h>

namespace hydrochrono::hydro {
// Forward declarations
class HydrostaticsComponent;
class RadiationComponent;
class ExcitationComponent;
class HydroSystem;
class ChronoHydroCoupler;
// Types GeneralizedForce and BodyForces are defined in hydro_types.h (included above)
}

using namespace chrono;
using namespace chrono::fea;

class ForceFunc6d;
class TestHydro;

class ComponentFunc : public ChFunction {
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
    ForceFunc6d* base_;  ///< Pointer to the full 6D force on the body.
    int index_;          ///< Index representing force degree of freedom on the body.
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
     * @param all_hydro_forces_user The TestHydro class where the total force on all bodies is calculated.
     */
    ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* all_hydro_forces_user);

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

    std::shared_ptr<ChBody> body_;                  ///< Pointer to the body this force is applied to.
    int b_num_;                                     ///< Body's index in the system. Currently 1-indexed.
    ComponentFunc forces_[6];                       ///< Forces for each degree of freedom.
    std::shared_ptr<ComponentFunc> force_ptrs_[6];  ///< Pointers to the forces.
    std::shared_ptr<ChForce> chrono_force_;         ///< Chrono force for the body.
    std::shared_ptr<ChForce> chrono_torque_;        ///< Chrono torque for the body.
    TestHydro* all_hydro_forces_;                   ///< Pointer to TestHydro for calculations.
};

class ChLoadAddedMass;

// Lightweight hydrodynamics profiling stats
struct HydroProfileStats {
    double hydrostatics_seconds = 0.0;
    double radiation_seconds    = 0.0;
    double waves_seconds        = 0.0;
    int hydrostatics_calls      = 0;
    int radiation_calls         = 0;
    int waves_calls             = 0;
};

// TODO: Rename TestHydro for clarity, perhaps to HydroForces?
// TODO: Split TestHydro class from its helper classes for clearer code structure.
class TestHydro {
  public:
    TestHydro() = delete;

    /**
     * @brief Main constructor for initializing the TestHydro class.
     *
     * Sets up vector of bodies, h5 file info, and hydro inputs. If no waves are given,
     * this constructor defaults to using NoWave.
     *
     * @param user_bodies List of pointers to bodies for the hydro forces.
     * @param h5_file_name Name of the h5 file where hydro data is stored.
     * @param waves WaveBase object. Defaults to NoWave if not provided.
     */
    TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
              std::string h5_file_name,
              std::shared_ptr<WaveBase> waves = std::make_shared<NoWave>());

    // Destructor (defined in .cpp to allow unique_ptr to incomplete type)
    ~TestHydro();

    // Deleted copy constructor and assignment operator for safety.
    TestHydro(const TestHydro& old) = delete;
    TestHydro& operator=(const TestHydro& rhs) = delete;

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
     * NOTE: Retained for backward compatibility. The main force path now uses HydroSystem
     * (via CoordinateFuncForBody). This method is not called during normal simulation.
     *
     * @return 6N dimensional force for 6 DOF and N bodies in system.
     */
    std::vector<double> ComputeForceHydrostatics();

    /**
     * @brief Legacy API: Computes radiation damping forces for a 6N dimensional system.
     *
     * NOTE: Retained for backward compatibility. The main force path now uses HydroSystem
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
     * NOTE: Retained for backward compatibility. The main force path now uses HydroSystem
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

    // Convolution mode selection
    enum class RadiationConvolutionMode {
        Baseline,
        TaperedDirect
    };

    /**
     * @brief Set the radiation convolution mode. Default is Baseline.
     */
    void SetRadiationConvolutionMode(RadiationConvolutionMode mode) {
        convolution_mode_ = mode;
        InvalidateRadiationComponent();  // Invalidate component to recreate with new settings
    }

        struct TaperedDirectOptions {
            // smoothing: "sg" (Savitzky–Golay) or "moving_average"
            std::string smoothing = "sg";
            int window_length = 5;                 // odd, >= 3
            
            // RIRF truncation
            double rirf_end_time = -1.0;           // end RIRF at this time (seconds), -1.0 = use full length
            
            // Simple taper control - sensible defaults for improved stability
            double taper_start_percent = 0.8;      // start taper at 80% (taper last 20%)
            double taper_end_percent = 1.0;        // end taper at 100% of total time series  
            double taper_final_amplitude = 0.0;    // final amplitude as fraction of original (0.0 = zero, 1.0 = no change)
            bool export_plot_csv = false;          // dump before/after CSV summaries (false by default)
        };

    /**
     * @brief Set options for TaperedDirect preprocessing.
     */
    void SetTaperedDirectOptions(const TaperedDirectOptions& opts) {
        tapered_opts_ = opts;
        InvalidateRadiationComponent();  // Invalidate component to recreate with new settings
    }

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
     * If the total force for the body and DOF was computed for the current timestep, it's retrieved.
     * Otherwise, the function calculates it. Note: Body index is 1-based here due to its origin from ForceFunc6d.
     *
     * @param b Body index (1-based due to source from ForceFunc6d).
     * @param i Degree of Freedom (DOF) index, ranging from [0,...5].
     *
     * @return Component of the force vector for body 'b' and DOF 'i'.
     */
    double CoordinateFuncForBody(int b, int i);

    // Hydrodynamics profiling accessors
    HydroProfileStats GetProfileStats() const { return profile_stats_; }

    // Compare mode: legacy debugging feature, retained for API compatibility.
    // No longer affects behavior since the main path now uses HydroSystem exclusively.
    void SetCompareMode(bool enable) { compare_mode_ = enable; }
    bool GetCompareMode() const { return compare_mode_; }

  private:
    // Class properties related to the body and hydrodynamics
    std::vector<std::shared_ptr<ChBody>> bodies_;
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

    // Cached SystemState: built once per time step and reused by all force computations
    hydrochrono::hydro::SystemState cached_state_;
    double cached_state_time_ = std::numeric_limits<double>::quiet_NaN();

    // Added mass related properties
    std::shared_ptr<ChLoadContainer> my_loadcontainer;
    std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;

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

    // Radiation damping force component
    std::unique_ptr<hydrochrono::hydro::RadiationComponent> radiation_component_;

    // Wave excitation force component
    std::unique_ptr<hydrochrono::hydro::ExcitationComponent> excitation_component_;

    // HydroSystem + ChronoHydroCoupler (persistent, constructed once)
    std::unique_ptr<hydrochrono::hydro::HydroSystem> hydro_system_;
    std::unique_ptr<hydrochrono::hydro::ChronoHydroCoupler> chrono_coupler_;

    // Legacy compare mode flag: retained for API compatibility; no longer used.
    bool compare_mode_ = false;

    // Convolution kernel preprocessing (optional)
    RadiationConvolutionMode convolution_mode_ = RadiationConvolutionMode::Baseline;
    bool rirf_processed_ready_ = false;
    std::vector<Eigen::Tensor<double, 3>> rirf_processed_; // per body [dof x col x step]
    TaperedDirectOptions tapered_opts_;
    std::string diagnostics_output_dir_;

    void EnsureProcessedRIRF();

    // Helper to ensure radiation component exists with current settings
    void EnsureRadiationComponent();
    // Helper to invalidate radiation component (requires full type definition)
    void InvalidateRadiationComponent();
    // Helper to ensure excitation component exists
    void EnsureExcitationComponent();

    // Factory: creates ExcitationComponent with current wave configuration.
    // Used by both EnsureExcitationComponent() and EnsureHydroSystemAndCoupler()
    // to ensure consistent construction.
    std::unique_ptr<hydrochrono::hydro::ExcitationComponent> CreateExcitationComponent() const;

    // Factory: creates HydrostaticsComponent with current equilibrium/stiffness data.
    // Used by both the TestHydro constructor and EnsureHydroSystemAndCoupler()
    // to ensure consistent construction.
    std::unique_ptr<hydrochrono::hydro::HydrostaticsComponent> CreateHydrostaticsComponent() const;

    // Factory: creates RadiationComponent with current BEM data and convolution settings.
    // Used by both EnsureRadiationComponent() and EnsureHydroSystemAndCoupler()
    // to ensure consistent construction.
    // Note: Each instance owns its own velocity history (they are independent).
    std::unique_ptr<hydrochrono::hydro::RadiationComponent> CreateRadiationComponent() const;

    // Internal helper: constructs hydro_system_ and chrono_coupler_ once.
    // Subsequent calls are no-ops. Called automatically by CoordinateFuncForBody().
    void EnsureHydroSystemAndCoupler();
};

#endif