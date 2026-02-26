---
layout: page
title: CMake Build Structure
parent_section: Developer Documentation
---

# HydroChrono CMake Build Structure

This page gives a brief overview of how `CMakeLists.txt` is organized, why it’s structured this way, and what to change when adding features. It reflects the current top‑level `CMakeLists.txt` in this repo.

This is internal developer documentation for contributors working on HydroChrono’s build/packaging.

---

## Overview

The top‑level file is split into clearly marked sections:

1. Prerequisites & Policies
2. Project Definition & Defaults
3. User Options
4. Find Dependencies
5. Core Library Target
6. Auxiliary Targets (GUI, Tests, Demos)
7. Installation & Packaging

Banner comments separate sections so it’s easy to navigate and reason about changes.

---

## 1) Prerequisites & Policies

- `cmake_minimum_required(3.18.2)` (CMP0091, modern target props)
- `CMP0091 NEW` to control MSVC runtime via `CMAKE_MSVC_RUNTIME_LIBRARY`
- Guard against in‑source builds

---

## 2) Project Definition & Defaults

- Project `HydroChrono` (C++)
- Default build type: Release (Debug/MinSizeRel/RelWithDebInfo supported)
- `cmake/` added to `CMAKE_MODULE_PATH` (local modules)
- Generates `version.h` under `build/hydroc/`

---

## 3) User Options (feature toggles)

- `HYDROCHRONO_ENABLE_TESTS` - Enable test targets
- `HYDROCHRONO_ENABLE_IRRLICHT` - Enable Irrlicht visualization
- `HYDROCHRONO_ENABLE_VSG` - Enable VSG visualization
- `HYDROCHRONO_ENABLE_DEMOS` - Enable demo executables
- `HYDROCHRONO_ENABLE_YAML_RUNNER` - Enable YAML-based CLI runner
- `HYDROCHRONO_ENABLE_MOORDYN` - Enable MoorDyn mooring coupling (requires `extern/MoorDyn` submodule)

Enable lean builds (e.g., CI) or developer variants.

---

## 4) Find Dependencies

### Chrono

- `find_package(Chrono CONFIG REQUIRED COMPONENTS Parsers)`
- Irrlicht and VSG targets are checked if the corresponding `HYDROCHRONO_ENABLE_*` option is ON
- HDF5 targets are pre-loaded before Chrono to avoid conflicts with `hdf5_tools` shim targets

### HDF5 & Eigen

- HDF5: first tries config mode (`find_package(HDF5 CONFIG QUIET)`), falls back to module mode (`find_package(HDF5 REQUIRED COMPONENTS CXX)`). Prefers static libraries when available.
- Eigen: auto-detected from Chrono's build configuration

### Platform notes

- MSVC: a few warning suppressions; Eigen aligned‑storage workaround

---

## 5) Core Library Target

Libraries:

- `HydroChrono` — core hydrodynamics (HDF5 I/O, YAML, radiation, waves, utilities)
- `HydroChronoGUI` — GUI helpers for Irrlicht and/or VSG visualization

Each target is configured directly with C++ standard, PIC, include dirs, and Chrono links.

Key links:

- `HydroChrono` → `Chrono::Chrono_core`, OpenMP, HDF5 (+ MoorDyn when enabled)
- `HydroChronoGUI` → `Chrono::Chrono_core` (+ `Chrono::Chrono_irrlicht` and/or `Chrono::Chrono_vsg` when enabled)


---

## 6) Auxiliary Targets (GUI, Tests, Demos)

### YAML‑Driven CLI

- Executable: `run_hydrochrono` (built on top of the `HydroChrono` library)
- Links: `HydroChrono`, `HydroChronoGUI`, `Chrono::Chrono_parsers`
- Minimal CLI wrapper: runner sources and small utils


>Why library + CLI split?
> - Clear separation: core in `HydroChrono`; thin CLI wrapper.
> - Reuse & testing: same core for tests/demos/GUI; tests link the lib directly.
> - Packaging & speed: multiple frontends, faster CLI relinks, keep extras out of core.

### Tests & Demos

- If enabled, tests are added via `add_subdirectory(tests/regression)` and `add_subdirectory(tests/unit)`
- A helper `configure_test_environment()` assembles PATH on Windows so Chrono, HDF5, Irrlicht, VSG, and MoorDyn DLLs are found when tests run from the build tree
- Demos can be included behind `HYDROCHRONO_ENABLE_DEMOS`

Runtime concerns that affect tests:
- Matching build types between Chrono and HydroChrono (Release vs Debug)
- DLL search order on Windows

---

## 7) Installation & Packaging

Two layers:

1) Dev kit (optional)
   - `HC_INSTALL_DEV_KIT` (OFF by default)
   - Installs headers, libs, and `HydroChronoConfig.cmake`/targets for downstream `find_package(HydroChrono)`

2) Runtime (default)
   - Flat layout via `cmake --install` + CPack ZIP
   - Installs `run_hydrochrono.exe`, required DLLs, curated tests
   - Copies Chrono visual assets into `data/`; includes `tests/RUN-TESTS.ps1`

CPack (ZIP): `cpack -C Release` yields a ready‑to‑share artifact.


---

## Tips for Modifying the Build

- Add new features behind options (e.g., `HYDROCHRONO_ENABLE_<FEATURE>`) with sensible defaults
- Update `configure_hydro_target()` instead of duplicating setup
- Prefer imported targets (e.g., `Chrono::Chrono_core`) over manual flags
- Avoid hardcoded paths; expose `CACHE` variables when needed
- Keep section banners to make diffs/reviews faster

---

## Related Docs

- [Build & Setup](build_instructions)
- [CMake Build Basics](cmake_build_basics)
- [Project Chrono Documentation](https://api.projectchrono.org/)

<p align="center">
  <img src="https://nrel.github.io/HydroChrono/assets/img/wave_animation2.gif" alt="Wave Energy" width="80%" />
</p>