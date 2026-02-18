// HydroChrono VSG Radiation Surface Visualization - Physics-Based Implementation
//
// Implements radiated wave visualization using linearized potential flow theory:
//
// Theory:
//   For a body oscillating in mode j with velocity ξ̇_j, the far-field radiated
//   wave elevation is approximately:
//
//     η_rad(r,θ,t) ≈ Σ_j  A_j · S_j(θ) · cos(kr - ωt + φ_j) · D(r)
//
//   Where:
//     - A_j = amplitude proportional to body velocity (ideally from Kochin functions)
//     - S_j(θ) = angular shape function (monopole for heave, dipole for surge/sway)
//     - k = wavenumber from dispersion relation ω² = gk·tanh(kh)
//     - D(r) = 1/√r cylindrical spreading decay
//
//   The amplitude envelope propagates at group velocity c_g, while wave crests
//   propagate at phase velocity c_p. For deep water: c_g = c_p/2.
//
// Limitations:
//   - Amplitude scaling is empirical (Kochin functions not available)
//   - Body geometry effects not captured beyond characteristic radius
//   - Linear superposition assumed (valid for small amplitudes)

#include "vsg_radiation_surface.h"

#include <algorithm>
#include <cmath>

namespace hydroc {
namespace gui {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kDispersionTolerance = 1e-10;
constexpr double kMinWavenumber = 1e-6;
}  // namespace

RadiationSurfaceViz::RadiationSurfaceViz() : params_() {
    UpdateWaveProperties();
}

RadiationSurfaceViz::RadiationSurfaceViz(const Params& p) : params_(p) {
    UpdateWaveProperties();
}

RadiationSurfaceViz::~RadiationSurfaceViz() = default;
RadiationSurfaceViz::RadiationSurfaceViz(RadiationSurfaceViz&&) noexcept = default;
RadiationSurfaceViz& RadiationSurfaceViz::operator=(RadiationSurfaceViz&&) noexcept = default;

void RadiationSurfaceViz::SetParams(const Params& p) {
    params_ = p;
    UpdateWaveProperties();
}

void RadiationSurfaceViz::Reset() {
    sources_.clear();
    current_time_ = 0.0;
}

// ----------------------------------------------------------------------------
// Wave property calculations
// ----------------------------------------------------------------------------

double RadiationSurfaceViz::SolveDispersion(double omega, double depth, double g, int max_iter) {
    // Solve ω² = gk·tanh(kh) for k using Newton-Raphson iteration.
    // Initial guess: deep water approximation k₀ = ω²/g.
    
    if (omega <= 0.0 || g <= 0.0) {
        return kMinWavenumber;
    }

    const double omega_sq = omega * omega;
    double k = omega_sq / g;  // Deep water initial guess

    // For very deep water (kh > ~3), tanh(kh) ≈ 1, so initial guess is good.
    // For finite depth, iterate.
    const double kh_deep = 3.0;
    if (k * depth < kh_deep) {
        // Newton-Raphson: f(k) = ω² - gk·tanh(kh) = 0
        // f'(k) = -g·tanh(kh) - gkh·sech²(kh)
        for (int i = 0; i < max_iter; ++i) {
            const double kh = k * depth;
            const double tanh_kh = std::tanh(kh);
            const double sech_kh = 1.0 / std::cosh(kh);
            
            const double f = omega_sq - g * k * tanh_kh;
            const double df = -g * tanh_kh - g * k * depth * sech_kh * sech_kh;
            
            if (std::abs(df) < kDispersionTolerance) {
                break;
            }
            
            const double dk = -f / df;
            k += dk;
            
            if (k < kMinWavenumber) {
                k = kMinWavenumber;
            }
            
            if (std::abs(dk) < kDispersionTolerance * k) {
                break;
            }
        }
    }

    return std::max(k, kMinWavenumber);
}

void RadiationSurfaceViz::UpdateWaveProperties() {
    const double T = std::max(params_.wave_period, 0.5);
    omega_ = kTwoPi / T;
    
    // Solve dispersion relation for wavenumber.
    k_ = SolveDispersion(omega_, params_.water_depth, params_.gravity);
    
    // Wavelength.
    lambda_ = kTwoPi / k_;
    
    // Phase speed: c_p = ω/k.
    c_phase_ = omega_ / k_;
    
    // Group speed: c_g = dω/dk = c_p · (1/2 + kh/sinh(2kh)) for finite depth.
    // For deep water (kh → ∞): c_g = c_p/2.
    // For shallow water (kh → 0): c_g = c_p.
    const double kh = k_ * params_.water_depth;
    double n;  // c_g/c_p ratio
    if (kh > 10.0) {
        // Deep water limit.
        n = 0.5;
    } else if (kh < 0.01) {
        // Shallow water limit.
        n = 1.0;
    } else {
        // General finite depth.
        n = 0.5 * (1.0 + 2.0 * kh / std::sinh(2.0 * kh));
    }
    c_group_ = n * c_phase_;
}

double RadiationSurfaceViz::GetWavenumber() const { return k_; }
double RadiationSurfaceViz::GetWavelength() const { return lambda_; }
double RadiationSurfaceViz::GetPhaseSpeed() const { return c_phase_; }
double RadiationSurfaceViz::GetGroupSpeed() const { return c_group_; }
double RadiationSurfaceViz::GetAngularFrequency() const { return omega_; }

// ----------------------------------------------------------------------------
// Motion history management
// ----------------------------------------------------------------------------

void RadiationSurfaceViz::PruneHistory(BodySource& source, double current_time) {
    const double cutoff = current_time - params_.history_duration;
    while (!source.history.empty() && source.history.front().time < cutoff) {
        source.history.pop_front();
    }
}

void RadiationSurfaceViz::SetSourceState(const std::string& body_id,
                                          const chrono::ChVector3d& pos,
                                          const chrono::ChVector3d& vel,
                                          const chrono::ChVector3d& ang_vel,
                                          double radius,
                                          double t) {
    current_time_ = t;

    // Get or create source for this body.
    BodySource& source = sources_[body_id];
    source.current_x = pos.x();
    source.current_y = pos.y();
    source.radius = radius;

    MotionSample sample;
    sample.time = t;
    sample.pos_x = pos.x();
    sample.pos_y = pos.y();
    sample.vel_x = vel.x();
    sample.vel_y = vel.y();
    sample.vel_z = vel.z();
    sample.ang_vel_x = ang_vel.x();  // Roll rate (rotation about X)
    sample.ang_vel_y = ang_vel.y();  // Pitch rate (rotation about Y)

    // Ensure monotonic time (reset history if time went backwards).
    if (!source.history.empty() && t <= source.history.back().time) {
        source.history.clear();
    }

    source.history.push_back(sample);
    PruneHistory(source, t);
}

RadiationSurfaceViz::MotionSample RadiationSurfaceViz::InterpolateMotion(
    const BodySource& source, double t_ret) const {
    
    MotionSample result{t_ret, 0, 0, 0, 0, 0, 0, 0};

    if (source.history.empty()) {
        return result;
    }

    // Before recorded history: return zero.
    if (t_ret <= source.history.front().time) {
        return result;
    }
    // After recorded history: return most recent.
    if (t_ret >= source.history.back().time) {
        return source.history.back();
    }

    // Binary search for bracketing samples.
    auto it = std::upper_bound(source.history.begin(), source.history.end(), t_ret,
        [](double t, const MotionSample& s) { return t < s.time; });

    if (it == source.history.begin() || it == source.history.end()) {
        return result;
    }

    const auto& s0 = *(it - 1);
    const auto& s1 = *it;

    const double dt = s1.time - s0.time;
    if (dt < 1e-10) {
        return s0;
    }
    const double alpha = (t_ret - s0.time) / dt;

    result.time = t_ret;
    result.pos_x = s0.pos_x + alpha * (s1.pos_x - s0.pos_x);
    result.pos_y = s0.pos_y + alpha * (s1.pos_y - s0.pos_y);
    result.vel_x = s0.vel_x + alpha * (s1.vel_x - s0.vel_x);
    result.vel_y = s0.vel_y + alpha * (s1.vel_y - s0.vel_y);
    result.vel_z = s0.vel_z + alpha * (s1.vel_z - s0.vel_z);
    result.ang_vel_x = s0.ang_vel_x + alpha * (s1.ang_vel_x - s0.ang_vel_x);
    result.ang_vel_y = s0.ang_vel_y + alpha * (s1.ang_vel_y - s0.ang_vel_y);

    return result;
}

// ----------------------------------------------------------------------------
// Wave elevation computation
// ----------------------------------------------------------------------------

double RadiationSurfaceViz::EvaluateBodyContribution(
    const BodySource& source, double x, double y, double t) const {
    
    if (source.history.empty()) {
        return 0.0;
    }

    // Vector from body to evaluation point.
    const double dx = x - source.current_x;
    const double dy = y - source.current_y;
    const double r_sq = dx * dx + dy * dy;
    const double r = std::sqrt(r_sq);
    const double theta = std::atan2(dy, dx);

    // Near-field cutoff based on body radius (not wavelength).
    const double r_min = source.radius * params_.near_field_radii;
    if (r < source.radius * 0.5) {
        return 0.0;  // Inside body.
    }

    // ------------------------------------
    // Retarded time for wave propagation (at phase velocity).
    // Using phase velocity gives more intuitive visual response where
    // wave crests appear to originate from the body position.
    // ------------------------------------
    const double t_ret = t - r / c_phase_;
    const MotionSample m = InterpolateMotion(source, t_ret);

    // ------------------------------------
    // Angular shape functions for each mode.
    // Based on far-field radiation patterns from potential flow:
    //   - Heave: monopole (uniform)          → S = 1
    //   - Surge: dipole aligned with x-axis  → S = cos(θ)
    //   - Sway:  dipole aligned with y-axis  → S = sin(θ)
    //   - Pitch: similar to surge (rotation about y-axis)
    //   - Roll:  similar to sway (rotation about x-axis)
    // ------------------------------------
    const double S_heave = 1.0;
    const double S_surge = std::cos(theta);
    const double S_sway = std::sin(theta);
    const double S_pitch = std::cos(theta);
    const double S_roll = std::sin(theta);

    // ------------------------------------
    // Effective velocity for each mode.
    // For rotational modes, convert angular velocity to equivalent linear
    // velocity at the waterplane using characteristic radius.
    // ------------------------------------
    const double r0 = source.radius;
    
    // Heave: vertical velocity directly creates waves.
    const double v_heave = m.vel_z;
    
    // Surge/Sway: horizontal velocities create directional waves.
    const double v_surge = m.vel_x;
    const double v_sway = m.vel_y;
    
    // Pitch/Roll: angular velocities create waves similar to surge/sway.
    // The effective "heave-like" velocity at the edge of the body is ω × r.
    const double v_pitch = m.ang_vel_y * r0;
    const double v_roll = m.ang_vel_x * r0;

    // ------------------------------------
    // Amplitude calculation.
    // In proper theory: A_j = (ω/g) · H_j(θ,ω) · ξ̇_j where H_j is Kochin function.
    // Without Kochin functions, we use an empirical relationship:
    //   A ∝ v · (ω/g) · √(r0)
    // The √(r0) factor accounts for larger bodies radiating more energy.
    //
    // Base amplification (empirical) ensures waves are visible for typical
    // body sizes and velocities. User can further adjust via visual_scale.
    // ------------------------------------
    constexpr double kBaseAmplification = 5.0;  // Empirical visibility boost
    const double amp_factor = kBaseAmplification * (omega_ / params_.gravity) * std::sqrt(r0) * params_.visual_scale;
    
    // Combine modal contributions with appropriate signs.
    // Negative heave velocity (moving down) creates positive wave crest ahead.
    const double A_heave = -v_heave * S_heave * amp_factor;
    const double A_surge = v_surge * S_surge * amp_factor * 0.7;  // Dipole slightly weaker
    const double A_sway = v_sway * S_sway * amp_factor * 0.7;
    const double A_pitch = -v_pitch * S_pitch * amp_factor * 0.5;
    const double A_roll = -v_roll * S_roll * amp_factor * 0.5;

    const double amplitude = A_heave + A_surge + A_sway + A_pitch + A_roll;

    // ------------------------------------
    // Spatial decay.
    // 1/√r for cylindrical spreading (2D wave propagation).
    // Additional exponential decay for visual effect (represents dissipation).
    // ------------------------------------
    const double decay_length = params_.decay_distance * lambda_;
    const double r_eff = std::max(r, r0);  // Avoid singularity at origin.
    
    // Cylindrical spreading: 1/√(kr) for r >> 1/k, but clamped for near-field.
    const double kr = k_ * r_eff;
    const double spreading = 1.0 / std::sqrt(std::max(kr, 0.5));
    
    // Exponential decay (visualization, represents viscous/breaking losses).
    const double exp_decay = std::exp(-r_eff / decay_length);
    
    // Near-field ramp: smooth transition from body edge to far-field.
    // Ramps from 0 at r=0.5*r0 to ~1 at r=r_min.
    const double ramp_start = r0 * 0.5;
    const double ramp_width = r_min - ramp_start;
    const double near_field_ramp = (ramp_width > 0.01) 
        ? std::min(1.0, std::max(0.0, (r - ramp_start) / ramp_width))
        : 1.0;

    // ------------------------------------
    // Wave structure from retarded time.
    // 
    // The oscillatory wave pattern emerges naturally from the retarded-time
    // velocity without needing an explicit cos(kr - ωt) term. When the body
    // oscillates at frequency ω_body:
    //
    //   v(t) = V₀ sin(ω_body · t)
    //   v(t_ret) = V₀ sin(ω_body · (t - r/c))
    //            = V₀ sin(ω_body·t - ω_body·r/c)
    //
    // With the dispersion relation c = g/ω (deep water), this becomes:
    //   v(t_ret) = V₀ sin(ω_body·t - k_body·r)
    //
    // So the wave pattern automatically matches the body's actual oscillation
    // frequency, not a prescribed wave_period parameter. This is essential for
    // decay tests where the body oscillates at its natural frequency.
    // ------------------------------------

    // ------------------------------------
    // Final wave elevation.
    // ------------------------------------
    return amplitude * spreading * exp_decay * near_field_ramp;
}

double RadiationSurfaceViz::EvaluateEta(double x, double y, double t) const {
    if (sources_.empty()) {
        return 0.0;
    }

    // Sum contributions from all bodies (linear superposition is valid
    // in linearized potential flow theory).
    double total_eta = 0.0;
    for (const auto& [body_id, source] : sources_) {
        total_eta += EvaluateBodyContribution(source, x, y, t);
    }

    return total_eta;
}

}  // namespace gui
}  // namespace hydroc
