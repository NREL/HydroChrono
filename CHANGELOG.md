# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
