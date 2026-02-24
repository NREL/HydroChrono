---
layout: page
title: Building HydroChrono
parent_section: Developer Documentation
---

# Build & Setup (Windows)

A PowerShell script (`build.ps1`) is provided to build HydroChrono from the command line. The steps to configure, build, and test with this script are outlined below. Alternatively, you can use the CMake GUI with Visual Studio; those instructions appear later on this page.

## Prerequisites

**Required:**

- Visual Studio 2022 (Desktop development with C++)
- CMake ≥ 3.18
- Project Chrono (built from source; enable the Parsers module at minimum)
- HDF5 (C++ libraries, any recent version with CMake config support)

**Auto-detected from Chrono** (you do not need to install these separately):

- Eigen (header-only linear algebra library, bundled or referenced by Chrono)
- Irrlicht / VSG (visualization backends, detected if Chrono was built with them)

**Optional:**

- Python 3.x (only needed for regression tests and building documentation)
- MoorDyn (mooring dynamics library, included as a git submodule -- see [MoorDyn setup](#moordyn-mooring-library) below)

## Quick build with build.ps1 (recommended)

1) Copy `build-config-example.json` to `build-config.json` and set paths:
```json
{
  "ChronoDir": "C:/path/to/chrono/build/cmake",
  "PythonRoot": "C:/path/to/python/env"
}
```

`ChronoDir` is the only required key. Eigen, HDF5, and visualization library paths are auto-detected from Chrono's build configuration. `PythonRoot` is optional and only needed if Python isn't on your PATH.

2) From the repo root, run (first build) - use `-Verbose` for a more detailed output:
```powershell
./build.ps1 -Verbose
```

3) Create a distributable ZIP (installs to `build/install` and zips it):
```powershell
./build.ps1 -Verbose -Package
```

### Switches

| Switch | Description |
|---|---|
| `-BuildType <type>` | `Release` (default), `Debug`, `RelWithDebInfo`, or `MinSizeRel` |
| `-Clean` | Remove existing `build/` before configuring |
| `-Verbose` | Show full CMake and MSBuild output |
| `-Package` | Run `cmake --install` and `cpack` to create a ZIP |
| `-NoIrrlicht` | Disable Irrlicht visualization (auto-enabled if Chrono has it) |
| `-NoVSG` | Disable VSG visualization (auto-enabled if Chrono has it) |
| `-MoorDyn` | Enable MoorDyn mooring coupling (requires submodule, see [below](#moordyn-mooring-library)) |
| `-NoDemos` | Skip building demo executables |
| `-NoTests` | Skip building test targets |
| `-ConfigPath <file>` | Use a different JSON config (default: `build-config.json`) |

Run `.\build.ps1 -Help` for the full reference.

Install tree layout:
```
install/
  bin/    # run_hydrochrono.exe + DLLs
  data/   # Chrono visual assets (skybox, colormaps)
  tests/  # public regression suite (run_hydrochrono)
```

## MoorDyn mooring library

[MoorDyn](https://github.com/FloatingArrayDesign/MoorDyn) is an optional mooring dynamics library included as a git submodule. To enable it:

1) Initialize the submodule (one-time setup, clones MoorDyn into `extern/MoorDyn`):
```powershell
git submodule update --init extern/MoorDyn
```

2) Build with the `-MoorDyn` flag:
```powershell
./build.ps1 -MoorDyn -Verbose
```

The submodule pins a specific MoorDyn version tested with HydroChrono. Running `git submodule update` after pulling new commits will keep it in sync.

## Alternative: CMake GUI + Visual Studio

1. Open CMake (GUI)
    - Source: HydroChrono repo root
    - Build: `<repo>/build`

2. Configure:
    - Set `Chrono_DIR` to `<chrono>/build/cmake` (required)
    - Eigen, Irrlicht, and VSG paths are auto-detected from Chrono's build
    - `HDF5_DIR` may need to be set if CMake doesn't find it automatically
    - Set `HYDROCHRONO_ENABLE_MOORDYN` to `ON` if using MoorDyn (requires submodule, see [above](#moordyn-mooring-library))
    - Ensure your `CMAKE_MSVC_RUNTIME_LIBRARY` matches the one used to build Chrono

3. Generate → Open in Visual Studio → Build `run_hydrochrono` (choose the matching config, e.g., Release)

4. Runtime assets
    - Copy required DLLs next to the exe, or use the `-Package` flow above to stage `install/bin`


## Run regression tests (from the source tree)

After a local build, you have two regression test suites - **Option A** uses CTest to test the C++ based HydroChrono models, and **Option B** uses Python to test the main HydroChrono app - `run_hydrochrono.exe` :

**Option A** - CTest:
```powershell
ctest -C Release -L regression --test-dir build
```

  - Common CTest options/examples:

```powershell
# Show failing test output
ctest -C Release -L regression --test-dir build --output-on-failure

# Extra verbose output
ctest -C Release -L regression --test-dir build -VV

# Only run specific tests
ctest -C Release -R f3of -L regression --test-dir build
ctest -C Release -R rm3  -L regression --test-dir build

# Run in parallel
ctest -C Release -j 6 -L regression --test-dir build
```

**Option B** — Python-based YAML tests for the main executable:
```powershell
cd .\tests\regression\run_hydrochrono
python .\run_tests.py --all --exe ..\..\..\build\bin\Release\run_hydrochrono.exe
# Run a subset:
# python .\run_tests.py --sphere-decay --exe ..\..\..\build\bin\Release\run_hydrochrono.exe
# Optional GUI: add --gui
```

Outputs:
- Results (HDF5) and plots are written per case under `tests\regression\run_hydrochrono\<case>\<test>\outputs\`.
- PASS/FAIL summary and RMSrel are printed to the console.

Note: For the packaged ZIP, use the included `tests\RUN-TESTS.ps1` or run `run_tests.py` with `--exe ..\..\bin\run_hydrochrono.exe`.

### Generate a PDF regression report

Create a consolidated PDF report of results and comparisons:
```powershell
python .\tests\regression\utilities\generate_report.py --build-dir .\build --pdf
```
This summarizes PASS/FAIL and embeds plots for the standard regression cases.

## Troubleshooting

- Chrono/HydroChrono build types must match (e.g., both Release, same version of Visual Studio, etc.)
- `yaml-cpp.dll` or other dlls missing → ensure they are next to the exe
- GUI skybox missing → ensure `data/skybox/` exists in the install ZIP (packaging now includes it)
- Python tests: install `numpy`, `h5py`, `PyYAML`, `matplotlib`, or run `tests/RUN-TESTS.ps1`

<p align="center">
  <img src="https://nrel.github.io/HydroChrono/assets/img/wave_animation2.gif" alt="Wave Energy" width="80%" />
</p>