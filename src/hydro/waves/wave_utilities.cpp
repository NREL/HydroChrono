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

#include <hydroc/logging.h>

bool is_in_deep_water(double wavenumber, double water_depth) {
    return (wavenumber * water_depth > 89.4);
}

Eigen::Vector3d GetWheelerStretchedPosition(const Eigen::Vector3d& position, double eta, double water_depth, double mwl) {
    // Position relative to mean water level
    double z_pos = position.z() - mwl;
    // Wheeler stretching
    double z_stretched = water_depth * (z_pos - eta) / (water_depth + eta);
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
    if (is_in_deep_water(wavenumber, water_depth)) {
        water_velocity[0] =
            omega * amplitude * std::exp(wavenumber * z_pos) * cos(wavenumber * x_pos - omega * time + phase);
        water_velocity[2] =
            omega * amplitude * std::exp(wavenumber * z_pos) * sin(wavenumber * x_pos - omega * time + phase);
    } else {
        water_velocity[0] = omega * amplitude * std::cosh(wavenumber * (z_pos + water_depth)) /
                            std::sinh(wavenumber * water_depth) * cos(wavenumber * x_pos - omega * time + phase);
        water_velocity[2] = omega * amplitude * std::sinh(wavenumber * (z_pos + water_depth)) /
                            std::sinh(wavenumber * water_depth) * sin(wavenumber * x_pos - omega * time + phase);
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
    if (is_in_deep_water(wavenumber, water_depth)) {
        water_acceleration[0] =
            omega * omega * amplitude * std::exp(wavenumber * z_pos) * sin(wavenumber * x_pos - omega * time + phase);
        water_acceleration[2] =
            -omega * omega * amplitude * std::exp(wavenumber * z_pos) * cos(wavenumber * x_pos - omega * time + phase);
    } else {
        water_acceleration[0] = omega * omega * amplitude * std::cosh(wavenumber * (z_pos + water_depth)) /
                                std::sinh(wavenumber * water_depth) * sin(wavenumber * x_pos - omega * time + phase);
        water_acceleration[2] = -omega * omega * amplitude * std::sinh(wavenumber * (z_pos + water_depth)) /
                                std::sinh(wavenumber * water_depth) * cos(wavenumber * x_pos - omega * time + phase);
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

std::vector<std::array<double, 3>> CreateFreeSurface3DPts(const Eigen::VectorXd& eta,
                                                          const Eigen::VectorXd& t_vec) {
    std::vector<std::array<double, 3>> surface(t_vec.size() * 2);

    for (int i = 0; i < t_vec.size(); ++i) {
        double t = -1 * t_vec[i];
        double z = eta[i];

        surface[2 * i]     = {t, -10.0, z};
        surface[2 * i + 1] = {t, 10.0, z};
    }

    return surface;
}

std::vector<std::array<size_t, 3>> CreateFreeSurfaceTriangles(size_t eta_size) {
    std::vector<std::array<size_t, 3>> triangles;

    for (size_t i = 0; i < eta_size / 2 - 1; ++i) {
        triangles.push_back({2 * i, 2 * i + 1, 2 * i + 3});
        triangles.push_back({2 * i, 2 * i + 3, 2 * i + 2});
    }

    return triangles;
}

void WriteFreeSurfaceMeshObj(const std::vector<std::array<double, 3>>& points,
                             const std::vector<std::array<size_t, 3>>& triangles,
                             const std::string& file_name) {
    std::ofstream out(file_name);
    if (!out) {
        std::cerr << "Failed to open " << file_name << std::endl;
        return;
    }

    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    out << "# Wavefront OBJ file exported by HydroChrono" << std::endl;
    out << "# File Created: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl << std::endl;

    out << "# Vertices: " << points.size() << std::endl << std::endl;
    out << std::fixed << std::setprecision(6);
    for (const auto& point : points) {
        out << "v ";
        out << std::setw(14) << point[0] << ' ';
        out << std::setw(14) << point[1] << ' ';
        out << std::setw(14) << point[2] << std::endl;
    }
    out << std::endl;

    out << "# Faces: " << triangles.size() << std::endl << std::endl;
    for (const auto& triangle : triangles) {
        out << "f ";
        out << std::setw(9) << triangle[0] + 1;
        out << std::setw(9) << triangle[1] + 1;
        out << std::setw(9) << triangle[2] + 1 << std::endl;
    }

    out.close();
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

