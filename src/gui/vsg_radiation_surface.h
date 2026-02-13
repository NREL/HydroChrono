// HydroChrono VSG Radiation Surface Visualization (Multi-Body, Retarded-Time)
//
// Implements a physics-based visualization of radiated waves from oscillating bodies.
// Uses linearized potential flow theory with the following features:
//   - Proper dispersion relation (deep water or finite depth)
//   - Group velocity for envelope propagation (energy transport)
//   - Phase velocity for wave crest propagation
//   - Angular radiation patterns: monopole (heave), dipole (surge/sway/pitch/roll)
//   - Cylindrical spreading (1/√r decay)
//
// NOTE: This is a visualization approximation. For accurate radiated wave fields,
// Kochin functions from BEM solvers would be required.
#pragma once

#include <chrono/core/ChVector3.h>
#include <deque>
#include <string>
#include <unordered_map>

namespace hydroc {
namespace gui {

/// Multi-body radiated wave visualization using linearized potential flow model.
/// Each body maintains its own motion history; contributions are linearly superposed.
class RadiationSurfaceViz {
  public:
    struct Params {
        double visual_scale = 1.0;       ///< Visual amplification multiplier
        
        /// Expected oscillation period for wave propagation (s).
        /// For wave-forced simulations: set to incident wave period.
        /// For decay tests: set to body's natural period T_n ≈ 2π√((M+A_∞)/K_hs).
        /// This affects the phase velocity used for wave propagation; mismatch
        /// between this and actual body frequency will distort the wave pattern.
        double wave_period = 8.0;
        
        double water_depth = 1000.0;     ///< Water depth (m); large values → deep water
        double gravity = 9.81;           ///< Gravitational acceleration (m/s²)
        double history_duration = 30.0;  ///< Motion history duration (s)
        double decay_distance = 4.0;     ///< Exponential decay distance in wavelengths
        double near_field_radii = 1.5;   ///< Near-field cutoff as multiple of body radius
    };

    RadiationSurfaceViz();
    explicit RadiationSurfaceViz(const Params& p);
    ~RadiationSurfaceViz();

    RadiationSurfaceViz(const RadiationSurfaceViz&) = delete;
    RadiationSurfaceViz& operator=(const RadiationSurfaceViz&) = delete;
    RadiationSurfaceViz(RadiationSurfaceViz&&) noexcept;
    RadiationSurfaceViz& operator=(RadiationSurfaceViz&&) noexcept;

    void SetParams(const Params& p);
    const Params& GetParams() const { return params_; }

    /// Update a source body's state (call for each body, each frame).
    /// @param body_id Unique identifier for the body (e.g., body name).
    /// @param pos Body center position in world coordinates.
    /// @param vel Body linear velocity in world coordinates.
    /// @param ang_vel Body angular velocity in world coordinates.
    /// @param radius Effective radiating radius (m) - affects near-field and amplitude.
    /// @param t Current simulation time.
    void SetSourceState(const std::string& body_id,
                        const chrono::ChVector3d& pos,
                        const chrono::ChVector3d& vel,
                        const chrono::ChVector3d& ang_vel,
                        double radius,
                        double t);

    /// Compute total radiated wave elevation from all bodies at position (x,y) and time t.
    /// @return Surface elevation η (m) - includes oscillatory wave structure.
    double EvaluateEta(double x, double y, double t) const;

    bool IsInitialized() const { return !sources_.empty(); }
    void Reset();

    /// Get number of active source bodies.
    size_t GetSourceCount() const { return sources_.size(); }

    /// Get computed wave properties (for debugging/display).
    double GetWavenumber() const;
    double GetWavelength() const;
    double GetPhaseSpeed() const;
    double GetGroupSpeed() const;
    double GetAngularFrequency() const;

  private:
    struct MotionSample {
        double time;
        double pos_x, pos_y;
        double vel_x, vel_y, vel_z;   // Linear velocity (m/s)
        double ang_vel_x, ang_vel_y;  // Roll/pitch rates (rad/s)
    };

    struct BodySource {
        std::deque<MotionSample> history;
        double current_x = 0.0;
        double current_y = 0.0;
        double radius = 5.0;
    };

    /// Interpolate motion at retarded time using stored history.
    MotionSample InterpolateMotion(const BodySource& source, double t_ret) const;

    /// Remove old samples beyond history_duration.
    void PruneHistory(BodySource& source, double current_time);

    /// Compute wave contribution from a single body at position (x,y) and time t.
    double EvaluateBodyContribution(const BodySource& source, double x, double y, double t) const;

    /// Recompute cached wave properties when params change.
    void UpdateWaveProperties();

    /// Solve dispersion relation ω² = gk·tanh(kh) for wavenumber k.
    static double SolveDispersion(double omega, double depth, double g, int max_iter = 50);

    Params params_;
    std::unordered_map<std::string, BodySource> sources_;
    double current_time_ = 0.0;

    // Cached wave properties (recomputed when params change).
    double omega_ = 0.0;       // Angular frequency (rad/s)
    double k_ = 0.0;           // Wavenumber (rad/m)
    double lambda_ = 0.0;      // Wavelength (m)
    double c_phase_ = 0.0;     // Phase speed (m/s)
    double c_group_ = 0.0;     // Group speed (m/s)
};

}  // namespace gui
}  // namespace hydroc
