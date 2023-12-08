---
layout: page
title: Building HydroChrono
parent_section: Developer Documentation
---

# Building HydroChrono (and demos)

## Introduction

This page provides instructions on how to build HydroChrono and its associated demos.

## Prerequisites

Before attempting to build HydroChrono, ensure you have installed and built all of the required [prerequisites](/developer_docs/prerequisites). The exact versions listed in the prerequisites are essential for proper functionality.

## Building HydroChrono

1. Clone this project into a directory and set up a build folder.
   - For instance, clone the project into a `HydroChrono` directory and set up a `HydroChrono_build` adjacent directory.

2. Using the CMake GUI, specify the location of the source files and built binaries for HydroChrono (`HydroChrono` and `HydroChrono_build` respectively). Configure and generate the solution with the following settings:
   
   - Set `Chrono_DIR` to the Chrono Build location (typically `../chrono_build/cmake`).
   - Set `HDF5_DIR` to your HDF5 build location, such as `../CMake-hdf5-1.10.8/CMake-hdf5-1.10.8/build/HDF5-1.10.8-win64/HDF5-1.10.8-win64/share/cmake`. Note that version 1.10.8 of HDF5 is best suited for Visual Studio 2019.
   - Enable the following options for additional features: `HYDROCHRONO_ENABLE_DEMOS`, `HYDROCHRONO_ENABLE_IRRLICHT`, and `HYDROCHRONO_ENABLE_TESTS`. For best results, enable all of these features. Note that the Irrlicht module requires Project Chrono to be built with the Irrlicht module enabled.
   - To build the docs: create a new entry in cmake-gui (Name: `Python3_ROOT_DIR`, Type: `PATH`) and set it to the path of your virtual Python environment where you have the `matplotlib`, `sphinx`, `sphinxcontrib-bibtex`, `breathe` and `h5py` packages installed (these are required to build the docs).

3. Navigate to the build folder and open the generated solution in Visual Studio (or click "Open Project" in the CMake GUI). Build the HydroChrono solution with your preferred configuration option (e.g. `RelWithDebInfo`). For comprehensive building and linking, use the `ALL_BUILD` project.

## Post-Build Steps

1. Copy the `chrono_build/bin/data` file from the Project Chrono build directory to `HydroChrono_build/data` to obtain optional shaders and logos.

2. Navigate to the `chrono_build/bin/RelWithDebInfo` directory and copy all `.dll` and `.pdb` files (excluding those for demos) to `HydroChrono_build/demos/RelWithDebInfo`. Files to copy include:

   - ChronoEngine.dll
   - ChronoEngine.pdb
   - ChronoEngine_irrlicht.dll
   - ChronoEngine_irrlicht.pdb
   - ChronoEngine_postprocess.dll
   - ChronoEngine_postprocess.pdb
   - Irrlicht.dll (note: no `.pdb` for Irrlicht)
   - Other `ChronoEngine_moduleName.dll` and `ChronoEngine_moduleName.pdb` files as relevant.

## Running Demos

1. To run the demos, navigate to `HydroChrono_build/demos/RelWithDebInfo`.

2. Demos require a command line argument indicating the location of input files. To specify this:
   
   - Run demo executables from the command line with an argument pointing to the absolute location of `<path>/HydroChrono/demos`. For example, executing `> ./sphere_decay.exe C:\\Users\\USERNAME\\HydroChrono\\demos` from the `HydroChrono_build/demos/RelWithDebInfo` directory would run the sphere heave decay test demo with the correct input file locations.
   - Alternatively, in Visual Studio 2019, you can set this command line argument for debug mode in a demo's properties. Navigate to the Solution Explorer, right-click a demo, select `properties > Configuration Properties > Debugging`, and set the Command Arguments to the appropriate path, like `C:\\Users\\USERNAME\\HydroChrono\\demos`.

3. Optionally, you can copy plot files into result directories and generate plots as needed.