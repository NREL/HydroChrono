/*********************************************************************
 * @file  wave_utilities.h
 * @brief Shared helper routines for wave models.
 *********************************************************************/

#ifndef HYDRO_WAVES_WAVE_UTILITIES_H
#define HYDRO_WAVES_WAVE_UTILITIES_H

#include <cstddef>
#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

bool is_in_deep_water(double wavenumber, double water_depth);

Eigen::Vector3d GetWheelerStretchedPosition(const Eigen::Vector3d& position,
                                            double                eta,
                                            double                water_depth,
                                            double                mwl);

double GetFreeSurfaceElevation(const Eigen::VectorXd& freqs_hz,
                               const Eigen::VectorXd& spectral_densities,
                               const Eigen::VectorXd& spectral_widths,
                               const Eigen::VectorXd& wave_phases,
                               const Eigen::VectorXd& wavenumbers,
                               const Eigen::Vector3d& position,
                               double time_value,
                               double water_depth);

std::vector<double> GetFreeSurfaceElevationTimeSeries(const Eigen::VectorXd& freqs_hz,
                                                      const Eigen::VectorXd& spectral_densities,
                                                      const Eigen::VectorXd& spectral_widths,
                                                      const Eigen::VectorXd& wave_phases,
                                                      const Eigen::VectorXd& wavenumbers,
                                                      const Eigen::Vector3d& position,
                                                      const Eigen::VectorXd& time_array,
                                                      double water_depth);

double GetEta(const Eigen::Vector3d& position,
              double time,
              double omega,
              double amplitude,
              double phase,
              double wavenumber);

double GetEtaIrregular(const Eigen::Vector3d& position,
                       double time,
                       const Eigen::VectorXd& freqs_hz,
                       const Eigen::VectorXd& spectral_densities,
                       const Eigen::VectorXd& spectral_widths,
                       const Eigen::VectorXd& wave_phases,
                       const Eigen::VectorXd& wavenumbers);

std::vector<double> GetEtaIrregularTimeSeries(const Eigen::Vector3d& position,
                                              const std::vector<double> time_index,
                                              const Eigen::VectorXd& freqs_hz,
                                              const Eigen::VectorXd& spectral_densities,
                                              const Eigen::VectorXd& spectral_widths,
                                              const Eigen::VectorXd& wave_phases,
                                              const Eigen::VectorXd& wavenumbers);

Eigen::Vector3d GetWaterVelocity(const Eigen::Vector3d& position,
                                 double time,
                                 double omega,
                                 double amplitude,
                                 double phase,
                                 double wavenumber,
                                 double water_depth,
                                 double mwl);

Eigen::Vector3d GetWaterAcceleration(const Eigen::Vector3d& position,
                                     double time,
                                     double omega,
                                     double amplitude,
                                     double phase,
                                     double wavenumber,
                                     double water_depth,
                                     double mwl);

Eigen::Vector3d GetWaterVelocityIrregular(const Eigen::Vector3d& position,
                                          double time,
                                          const Eigen::VectorXd& freqs_hz,
                                          const Eigen::VectorXd& spectral_densities,
                                          const Eigen::VectorXd& spectral_widths,
                                          const Eigen::VectorXd& wave_phases,
                                          const Eigen::VectorXd& wavenumbers,
                                          double water_depth,
                                          double mwl);

Eigen::Vector3d GetWaterAccelerationIrregular(const Eigen::Vector3d& position,
                                              double time,
                                              const Eigen::VectorXd& freqs_hz,
                                              const Eigen::VectorXd& spectral_densities,
                                              const Eigen::VectorXd& spectral_widths,
                                              const Eigen::VectorXd& wave_phases,
                                              const Eigen::VectorXd& wavenumbers,
                                              double water_depth,
                                              double mwl);

double ComputeWaveNumber(double omega,
                         double water_depth,
                         double g,
                         double tolerance = 1e-6,
                         int max_iterations = 100);

Eigen::VectorXd ComputeWaveNumbers(const Eigen::VectorXd& omegas,
                                   double water_depth,
                                   double g,
                                   double tolerance = 1e-6,
                                   int max_iterations = 100);

Eigen::VectorXd PiersonMoskowitzSpectrumHz(Eigen::VectorXd& f, double Hs, double Tp);

Eigen::VectorXd JONSWAPSpectrumHz(Eigen::VectorXd& f,
                                  double Hs,
                                  double Tp,
                                  double gamma       = 3.3,
                                  bool is_normalized = false);

std::vector<std::array<double, 3>> CreateFreeSurface3DPts(const Eigen::VectorXd& eta,
                                                          const Eigen::VectorXd& t_vec);

std::vector<std::array<size_t, 3>> CreateFreeSurfaceTriangles(size_t eta_size);

void WriteFreeSurfaceMeshObj(const std::vector<std::array<double, 3>>& points,
                             const std::vector<std::array<size_t, 3>>& triangles,
                             const std::string& file_name);

Eigen::VectorXd GetWidthArray(const Eigen::VectorXd& input_array);

#endif  // HYDRO_WAVES_WAVE_UTILITIES_H

