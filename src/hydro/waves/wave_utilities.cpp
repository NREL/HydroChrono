/*********************************************************************
 * @file  wave_utilities.cpp
 * @brief Helper routines shared by wave models.
 *********************************************************************/

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "wave_utilities.h"

#include <cmath>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#include <hydroc/logging.h>

bool is_in_deep_water(double wavenumber, double water_depth) {
    // Keep behavior consistent with ComputeWaveNumber(), which treats non-positive depth, very large
    // depth, or infinite depth as "effectively infinite".
    constexpr double EFFECTIVE_INFINITE_DEPTH = 1000.0;
    if (water_depth <= 0.0 || water_depth > EFFECTIVE_INFINITE_DEPTH || std::isinf(water_depth)) {
        return true;
    }

    // Numerical threshold for using deep-water asymptotic form (tanh(kh) ≈ 1).
    return (std::abs(wavenumber) * water_depth > 89.4);
}

Eigen::Vector3d GetWheelerStretchedPosition(const Eigen::Vector3d& position, double eta, double water_depth, double mwl) {
    // Wheeler stretching assumes finite depth. If depth is being used as a sentinel for "infinite",
    // fall back to no stretching to avoid undefined behavior.
    if (water_depth <= 0.0 || std::isinf(water_depth)) {
        return position;
    }

    // Position relative to mean water level
    double z_pos = position.z() - mwl;
    // Wheeler stretching
    double denom = water_depth + eta;
    const double scale = std::max(1.0, std::abs(water_depth));
    const double eps   = std::max(1e-12 * scale, 100.0 * std::numeric_limits<double>::epsilon() * scale);
    if (std::abs(denom) < eps) {
        // Pathological case (e.g., eta ≈ -water_depth): avoid division-by-zero blow-up.
        return position;
    }
    double z_stretched = water_depth * (z_pos - eta) / denom;
    return Eigen::Vector3d(position.x(), position.y(), z_stretched + mwl);
}

double GetEta(const Eigen::Vector3d& position,
              double time,
              double omega,
              double amplitude,
              double phase,
              double wavenumber) {
    double x_pos = position.x();
    return amplitude * cos(wavenumber * x_pos - omega * time + phase);
}

double GetEtaGradientX(const Eigen::Vector3d& position,
                       double time,
                       double omega,
                       double amplitude,
                       double phase,
                       double wavenumber) {
    // Analytic derivative of η = A·cos(k·x − ω·t + φ):
    //   ∂η/∂x = −A·k·sin(k·x − ω·t + φ)
    double x_pos = position.x();
    return -amplitude * wavenumber * sin(wavenumber * x_pos - omega * time + phase);
}

double GetEtaIrregular(const Eigen::Vector3d& position,
                       double time,
                       const Eigen::VectorXd& freqs_hz,
                       const Eigen::VectorXd& spectral_densities,
                       const Eigen::VectorXd& spectral_widths,
                       const Eigen::VectorXd& wave_phases,
                       const Eigen::VectorXd& wavenumbers) {
    double eta = 0.0;
    for (int i = 0; i < freqs_hz.size(); ++i) {
        double amplitude = std::sqrt(2 * spectral_densities[i] * spectral_widths[i]);
        double omega     = 2 * M_PI * freqs_hz[i];
        eta += GetEta(position, time, omega, amplitude, wave_phases[i], wavenumbers[i]);
    }
    return eta;
}

double GetEtaGradientXIrregular(const Eigen::Vector3d& position,
                                double time,
                                const Eigen::VectorXd& freqs_hz,
                                const Eigen::VectorXd& spectral_densities,
                                const Eigen::VectorXd& spectral_widths,
                                const Eigen::VectorXd& wave_phases,
                                const Eigen::VectorXd& wavenumbers) {
    // Sum of component gradients: ∂η/∂x = Σ[−Aᵢ·kᵢ·sin(kᵢ·x − ωᵢ·t + φᵢ)]
    double deta_dx = 0.0;
    for (Eigen::Index i = 0; i < freqs_hz.size(); ++i) {
        double amplitude = std::sqrt(2.0 * spectral_densities[i] * spectral_widths[i]);
        double omega     = 2.0 * M_PI * freqs_hz[i];
        deta_dx += GetEtaGradientX(position, time, omega, amplitude, wave_phases[i], wavenumbers[i]);
    }
    return deta_dx;
}

std::vector<double> GetEtaIrregularTimeSeries(const Eigen::Vector3d& position,
                                              const std::vector<double> time_index,
                                              const Eigen::VectorXd& freqs_hz,
                                              const Eigen::VectorXd& spectral_densities,
                                              const Eigen::VectorXd& spectral_widths,
                                              const Eigen::VectorXd& wave_phases,
                                              const Eigen::VectorXd& wavenumbers) {
    std::vector<double> eta(time_index.size(), 0.0);
    for (size_t j = 0; j < time_index.size(); ++j) {
        eta[j] = GetEtaIrregular(position,
                                 time_index[j],
                                 freqs_hz,
                                 spectral_densities,
                                 spectral_widths,
                                 wave_phases,
                                 wavenumbers);
    }
    return eta;
}

Eigen::Vector3d GetWaterVelocity(const Eigen::Vector3d& position,
                                 double time,
                                 double omega,
                                 double amplitude,
                                 double phase,
                                 double wavenumber,
                                 double water_depth,
                                 double mwl) {
    double x_pos = position.x();
    double z_pos = position.z() - mwl;

    Eigen::Vector3d water_velocity(0.0, 0.0, 0.0);
    const double k_mag = std::abs(wavenumber);
    if (!std::isfinite(k_mag) || k_mag == 0.0 || !std::isfinite(omega) || omega == 0.0 || !std::isfinite(amplitude) || amplitude == 0.0) {
        return water_velocity;
    }

    const double phase_arg = wavenumber * x_pos - omega * time + phase;

    if (is_in_deep_water(wavenumber, water_depth)) {
        water_velocity[0] =
            omega * amplitude * std::exp(k_mag * z_pos) * cos(phase_arg);
        water_velocity[2] =
            omega * amplitude * std::exp(k_mag * z_pos) * sin(phase_arg);
    } else {
        const double kh = k_mag * water_depth;
        if (std::abs(kh) < 1e-8) {
            // Small-kh safe form (avoids sinh(kh)→0 leading to inf*0 -> NaN).
            // Uses the limiting behavior: cosh(k(z+h))/sinh(kh) ~ 1/(kh) and sinh(k(z+h))/sinh(kh) ~ (z+h)/h.
            const double omega_over_k = omega / k_mag;
            water_velocity[0]         = omega_over_k * amplitude / water_depth * cos(phase_arg);
            water_velocity[2]         = omega * amplitude * ((z_pos + water_depth) / water_depth) * sin(phase_arg);
        } else {
            const double denom = std::sinh(kh);
            water_velocity[0] = omega * amplitude * std::cosh(k_mag * (z_pos + water_depth)) / denom * cos(phase_arg);
            water_velocity[2] = omega * amplitude * std::sinh(k_mag * (z_pos + water_depth)) / denom * sin(phase_arg);
        }
    }
    return water_velocity;
}

Eigen::Vector3d GetWaterAcceleration(const Eigen::Vector3d& position,
                                     double time,
                                     double omega,
                                     double amplitude,
                                     double phase,
                                     double wavenumber,
                                     double water_depth,
                                     double mwl) {
    double x_pos = position.x();
    double z_pos = position.z() - mwl;

    Eigen::Vector3d water_acceleration(0.0, 0.0, 0.0);
    const double k_mag = std::abs(wavenumber);
    if (!std::isfinite(k_mag) || k_mag == 0.0 || !std::isfinite(omega) || omega == 0.0 || !std::isfinite(amplitude) || amplitude == 0.0) {
        return water_acceleration;
    }

    const double phase_arg = wavenumber * x_pos - omega * time + phase;

    if (is_in_deep_water(wavenumber, water_depth)) {
        water_acceleration[0] =
            omega * omega * amplitude * std::exp(k_mag * z_pos) * sin(phase_arg);
        water_acceleration[2] =
            -omega * omega * amplitude * std::exp(k_mag * z_pos) * cos(phase_arg);
    } else {
        const double kh = k_mag * water_depth;
        if (std::abs(kh) < 1e-8) {
            // Small-kh safe form (avoids sinh(kh)→0 leading to inf*0 -> NaN).
            const double omega_over_k = omega / k_mag;
            water_acceleration[0]     = omega * omega_over_k * amplitude / water_depth * sin(phase_arg);
            water_acceleration[2]     = -omega * omega * amplitude * ((z_pos + water_depth) / water_depth) * cos(phase_arg);
        } else {
            const double denom = std::sinh(kh);
            water_acceleration[0] = omega * omega * amplitude * std::cosh(k_mag * (z_pos + water_depth)) / denom * sin(phase_arg);
            water_acceleration[2] = -omega * omega * amplitude * std::sinh(k_mag * (z_pos + water_depth)) / denom * cos(phase_arg);
        }
    }
    return water_acceleration;
}

Eigen::Vector3d GetWaterVelocityIrregular(const Eigen::Vector3d& position,
                                          double time,
                                          const Eigen::VectorXd& freqs_hz,
                                          const Eigen::VectorXd& spectral_densities,
                                          const Eigen::VectorXd& spectral_widths,
                                          const Eigen::VectorXd& wave_phases,
                                          const Eigen::VectorXd& wavenumbers,
                                          double water_depth,
                                          double mwl) {
    Eigen::Vector3d water_velocity(0.0, 0.0, 0.0);
    for (int i = 0; i < freqs_hz.size(); ++i) {
        double amplitude = std::sqrt(2 * spectral_densities[i] * spectral_widths[i]);
        double omega     = 2 * M_PI * freqs_hz[i];
        water_velocity += GetWaterVelocity(position, time, omega, amplitude, wave_phases[i], wavenumbers[i], water_depth, mwl);
    }
    return water_velocity;
}

Eigen::Vector3d GetWaterAccelerationIrregular(const Eigen::Vector3d& position,
                                              double time,
                                              const Eigen::VectorXd& freqs_hz,
                                              const Eigen::VectorXd& spectral_densities,
                                              const Eigen::VectorXd& spectral_widths,
                                              const Eigen::VectorXd& wave_phases,
                                              const Eigen::VectorXd& wavenumbers,
                                              double water_depth,
                                              double mwl) {
    Eigen::Vector3d water_acceleration(0.0, 0.0, 0.0);
    for (int i = 0; i < freqs_hz.size(); ++i) {
        double amplitude = std::sqrt(2 * spectral_densities[i] * spectral_widths[i]);
        double omega     = 2 * M_PI * freqs_hz[i];
        water_acceleration +=
            GetWaterAcceleration(position, time, omega, amplitude, wave_phases[i], wavenumbers[i], water_depth, mwl);
    }
    return water_acceleration;
}

double ComputeWaveNumber(double omega,
                         double water_depth,
                         double g,
                         double tolerance,
                         int max_iterations) {
    constexpr double DEEP_WATER_THRESHOLD = 1000.0;

    LOG_DEBUG("Computing wave number for omega=" << omega << ", water_depth=" << water_depth << ", g=" << g);

    if (omega <= 0.0) {
        LOG_ERROR("Invalid angular frequency: " << omega);
        throw std::runtime_error("Angular frequency must be positive.");
    }
    if (water_depth < 0.0) {
        LOG_ERROR("Invalid water depth: " << water_depth);
        throw std::runtime_error("Water depth cannot be negative.");
    }
    if (g <= 0.0) {
        LOG_ERROR("Invalid gravity: " << g);
        throw std::runtime_error("Gravity must be positive.");
    }
    if (tolerance <= 0.0) {
        LOG_ERROR("Invalid tolerance: " << tolerance);
        throw std::runtime_error("Tolerance must be positive.");
    }
    if (max_iterations <= 0) {
        LOG_ERROR("Invalid max iterations: " << max_iterations);
        throw std::runtime_error("Maximum iterations must be positive.");
    }

    if (water_depth == 0.0 || water_depth > DEEP_WATER_THRESHOLD || std::isinf(water_depth)) {
        LOG_DEBUG("Using deep water approximation");
        return omega * omega / g;
    }

    double k          = omega * omega / g;
    int    iterations = 0;
    double error      = 1.0;
    while (error > tolerance && iterations < max_iterations) {
        double tanh_kh = std::tanh(k * water_depth);
        double f       = omega * omega - g * k * tanh_kh;
        double df      = -2.0 * g * tanh_kh - g * k * water_depth * (1.0 - tanh_kh * tanh_kh);

        if (std::abs(df) < tolerance) {
            LOG_ERROR("Numerical instability: derivative too close to zero (df=" << df << ")");
            throw std::runtime_error("Numerical instability: derivative too close to zero.");
        }

        double delta_k = f / df;
        k -= delta_k;
        error = std::abs(delta_k);
        iterations++;
    }

    if (iterations >= max_iterations) {
        LOG_ERROR("Failed to converge after " << iterations << " iterations (final error=" << error << ")");
        throw std::runtime_error("Failed to converge within maximum iterations.");
    }

    LOG_DEBUG("Wave number computation complete: k=" << k << " (" << iterations << " iterations)");
    return k;
}

Eigen::VectorXd ComputeWaveNumbers(const Eigen::VectorXd& omegas,
                                   double water_depth,
                                   double g,
                                   double tolerance,
                                   int max_iterations) {
    Eigen::VectorXd wavenumbers(omegas.size());
    for (int i = 0; i < omegas.size(); ++i) {
        wavenumbers[i] = ComputeWaveNumber(omegas[i], water_depth, g, tolerance, max_iterations);
    }
    return wavenumbers;
}

Eigen::VectorXd PiersonMoskowitzSpectrumHz(Eigen::VectorXd& f, double Hs, double Tp) {
    std::sort(f.begin(), f.end());

    Eigen::VectorXd spectral_densities(f.size());

    for (int i = 0; i < f.size(); ++i) {
        spectral_densities[i] = 1.25 * std::pow(1 / Tp, 4) * std::pow(Hs / 2, 2) * std::pow(f[i], -5) *
                                std::exp(-1.25 * std::pow(1 / Tp, 4) * std::pow(f[i], -4));
    }

    return spectral_densities;
}

Eigen::VectorXd JONSWAPSpectrumHz(Eigen::VectorXd& f,
                                  double Hs,
                                  double Tp,
                                  double gamma,
                                  bool is_normalized) {
    auto spectral_densities = PiersonMoskowitzSpectrumHz(f, Hs, Tp);

    double normalization_factor = (1 - 0.287 * std::log(gamma));

    for (int i = 0; i < spectral_densities.size(); ++i) {
        double sigma = (f[i] <= 1.0 / Tp) ? 0.07 : 0.09;
        spectral_densities[i] *=
            std::pow(gamma, std::exp(-(1.0 / (2.0 * std::pow(sigma, 2))) * std::pow(f[i] * Tp - 1.0, 2)));
        if (is_normalized) {
            spectral_densities[i] *= normalization_factor;
        }
    }
    return spectral_densities;
}

Eigen::VectorXd GetWidthArray(const Eigen::VectorXd& input_array) {
    Eigen::VectorXd width_array(input_array.size());
    for (int ii = 0; ii < width_array.size(); ii++) {
        width_array[ii] = 0.0;
        if (ii < input_array.size() - 1) {
            width_array[ii] += 0.5 * std::abs(input_array[ii + 1] - input_array[ii]);
        }
        if (ii > 0) {
            width_array[ii] += 0.5 * std::abs(input_array[ii] - input_array[ii - 1]);
        }
    }
    return width_array;
}

