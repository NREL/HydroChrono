# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] — 2026-02-19

### Added

- **State-space radiation damping** — alternative to convolution-based radiation that approximates RIRF kernels as a sum of exponential and oscillatory modes, reducing per-timestep cost from O(N) to O(1). New classes: `RadiationStateSpaceModel`, `RadiationStateSpaceFitter`, `RadiationStateSpaceComponent`
- **`RadiationMethod` enum and `StateSpaceOptions` configuration struct** for selecting between convolution and state-space radiation at runtime
- **YAML configuration keys** — `radiation_method`, `ss_max_order`, `ss_r2_threshold`, `output_kernel_fit`
- **Kernel-fit diagnostics** — R² values and mode counts per DOF pair, optionally written to HDF5 output via `output_kernel_fit: true`
- **Unit tests** for `RadiationStateSpaceModel` and `RadiationStateSpaceFitter`
- **Regression tests** comparing state-space vs convolution for sphere decay, OSWEC decay, sphere irregular waves, and OSWEC irregular waves
- **Frequency-domain excitation transfer function** for irregular waves — replaces time-domain convolution with a pre-computed transfer function, evaluating in O(N_freq) per DOF per timestep instead of O(N_irf × N_freq)

### Changed

- **Irregular wave model decoupled from simulation parameters** — removed `simulation_dt_` and `simulation_duration_` from `IrregularWaveParams`, eliminating the dependency that prevented adaptive time integration. `nfrequencies_` is now a required parameter.
- **Removed redundant `num_bodies` from wave classes** — ownership moved to `WaveBase`, set automatically in `HydroSystem::AddWaves()`, simplifying the wave creation API and eliminating a source of configuration errors
- Wave class code quality cleanup

### Removed

- Stale `include/hydroc/hydro_forces.h` (replaced by `hydro_system.h`)
- Duplicate demo `demos/sphere/sphere_irreg_waves.cpp`

### Fixed

- **State-space fitter: order-1 fits were rejected** — the fitting loop started at order 2 and `FitFromSVD` explicitly rejected order < 2, so single-exponential kernels (rank-1 Hankel matrix) always returned invalid results

## [0.5.0] — 2026-02-16

### Added

- **VSG free-surface visualization** — animated water-surface mesh driven by the wave model, with per-vertex normals and configurable grid resolution (`vsg_water_surface`)
- **Radiated wave overlay** — visualization of body-generated wave patterns layered on the free surface (`vsg_radiation_surface`)
- **PBR ocean materials, scene lighting, and GUI stats component** — new modules `vsg_materials`, `vsg_lighting`, `vsg_gui_component`, `vsg_config`
- **Wave model API extensions** — `GetElevationGradientXY()` and `GetElevationForVisualization()` added to `WaveBase`, `RegularWave`, and `IrregularWave`
- **Added-mass determinism unit test** — `test_added_mass_determinism` verifies bit-identical added-mass assembly across independent trials and sweeps all Chrono solver types
- **Regression tests produce HDF5 output** validated by external Python comparison scripts; tests now run fully headless (all GUI/visualization code stripped from test builds)
- **CPack packaging** — project DLLs (`HydroChrono`, `HydroChronoGUI`) now included in the installer ZIP
- **OpenMP runtime** (`vcomp140.dll`) explicitly installed for MSVC packages
- **Chrono DLL collection** changed to glob all DLLs from Chrono's bin directory, capturing transitive VSG/yaml-cpp/draco dependencies

### Changed

- **Switched added-mass implementation to Chrono's built-in `ChLoadHydrodynamics`** — replaces HydroChrono's legacy `ChLoadAddedMass` (`ChLoadCustomMultiple`-based) as the default. The legacy implementation is retained in the codebase and can be re-activated by uncommenting `#define HYDROCHRONO_USE_LEGACY_ADDED_MASS` in `hydro_system.h`.
- Default solver changed from GMRES to SPARSE_QR for most regression tests (faster, deterministic)
- Irregular wave surface evaluation performance improved (parallelism and caching)
- Regression test report generator now correctly parses CTest logs for PASS/FAIL status instead of assuming PASS when plots exist

### Fixed

- **YAML runner: `LoadSolverData` never called** — YAML structure mismatch prevented solver data from being loaded
- **YAML runner: mesh file paths broken** — `m_script_directory` was empty, breaking relative `model_file:` paths
- **Regression report image paths** — report generator now computes correct relative paths to plot images instead of using hardcoded paths
- OSWEC solver switched from SPARSE_QR to GMRES to prevent divergence (see Known Issues)
- RM3 decay test fixes and cleanup
- Sphere irregular wave test default arguments corrected
- **Regular wave excitation phase indexing for multi-body systems** — `RegularWave::GetForceAtTime()` used `excitation_force_phase_[rowEx]` instead of `excitation_force_phase_[body_offset + rowEx]`, causing body 2+ to use body 1's excitation phase. This affected the RM3 plate (and any second+ body) heave response in regular waves. Single-body models (e.g., sphere) were unaffected.

### Known Issues

- **OSWEC tests use GMRES solver.** OSWEC demos and tests use GMRES rather than SPARSE_QR. SPARSE_QR may work for OSWEC but has not yet been validated.
- **PSOR solver cannot handle added-mass assembly.** The added-mass determinism unit test solver sweep reports that PSOR errors out because it cannot handle stiffness/damping matrices. All other swept solvers (SPARSE_QR, SPARSE_LU, MINRES, GMRES, BICGSTAB, BARZILAIBORWEIN, APGD) pass assembly.

## [0.4.0] — 2025

- YAML-driven CLI (`run_hydrochrono`) for running simulations from text-based configuration files
- Cummins-equation time-domain solver with BEM hydrodynamic coefficients (BEMIO HDF5 format)
- Multibody system support via Project Chrono (bodies, joints, actuators)
- Irrlicht run-time visualization (optional)
- HDF5 output for post-processing and validation
- Regression test suite (IEA sphere, OSWEC, RM3, F3OF / DeepCWind)
- Windows installer (ZIP) via CPack
