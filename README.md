# HydroChrono

**HydroChrono** is an emerging hydrodynamics simulation tool designed to model complex ocean systems. Seamlessly integrated with the [Project Chrono](https://projectchrono.org/) physics engine, it offers a powerful C++ API for a wide range of simulations.

## Capabilities

HydroChrono extends the versatility of Project Chrono to hydrodynamic applications, enabling you to simulate multibody wave energy converters (WECs) and floating offshore wind turbines (FOWTs) with Chrono::FEA for flexible bodies.

## Get Started

Visit the [HydroChrono website](https://nrel.github.io/HydroChrono/) for information about the underlying theory, build instructions, and usage details. Currently, HydroChrono is available as a source build, offering customization and optimization opportunities. A binary distribution with Python wrapper is in development and coming soon.

## Limitations

- Currently only supports first-order linear potential flow forces.
- Hydrodynamic forces are currently limited to rigid bodies.

## Vision and Future Goals

- Expand hydrodynamics to incorporate nonlinear and 2nd order forces.
- Integrate advanced simulation features using Project Chrono's DEM and FEA modules.
- Support seamless transitioning from potential flow to CFD and SPH for detailed FSI analysis.
- Develop a Python API to broaden accessibility and ease of use.
- Foster and support open-source collaboration in the research community.

## Build Instructions

### Prerequisites

- CMake 3.16 or higher
- A C++ compiler (Visual Studio 2019 or higher, or GCC through MinGW/MSYS2)
- Project Chrono (built from source, tested with v9.0.0 and v9.0.1)
- HDF5 1.10.8 or higher
- Python 3.8 or higher (with numpy and matplotlib)

### Building from Source

1. Clone the repository:
```powershell
git clone https://github.com/NREL/HydroChrono.git
cd HydroChrono
```

2. Create a build directory and navigate to it:
```powershell
mkdir build
cd build
```

3. Configure and build the project:
```powershell
cmake .. -DChrono_DIR="<path_to_chrono_build>/cmake" -DHDF5_DIR="<path_to_hdf5_cmake>" -DPython3_ROOT_DIR="<path_to_python>"
cmake --build . --config Release
```

4. Run the tests:
```powershell
ctest -C Release --output-on-failure
```

### Clean Build

To perform a clean rebuild of the project:
```powershell
# From the project root
Remove-Item -Recurse -Force build
mkdir build
cd build
```

### Post-Build Steps

Copy the following DLL files from your Chrono build directory to your build directory's `demos/Release` folder:
- ChronoEngine.dll
- ChronoEngine_irrlicht.dll (if using Irrlicht)
- Irrlicht.dll (if using Irrlicht)

For detailed build instructions, including Visual Studio setup and running demos, see the [developer documentation](https://nrel.github.io/HydroChrono/developer_docs/build_instructions.html).
