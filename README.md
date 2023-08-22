# Building HydroChrono with Project Chrono and HDF5 files

## Prerequisites 
1. Chrono installed and built in RelWithDebInfo mode v8.0.0 (built with Irrlicht module) and Chrono dependencies (see [Chrono Install/Build Guide for details](https://api.projectchrono.org/tutorial_install_chrono.html) ). This includes items like
	* C++ Compiler (Visual Studio 2019 or newer required, 2019 recommended)
	* Eigen3 (3.3.0 or newer required, 3.4.0 recommended)
	* Irrlicht Engine (1.8.4)
	* CMake
	* GIT client

2. Install H5Cpp header file from HDF5Group Note: Version 1.10.8 recommended, other versions may not work as intended. [Download Link](https://portal.hdfgroup.org/display/support/Downloads) [Instructions](https://portal.hdfgroup.org/display/support/Building+HDF5+with+CMake#BuildingHDF5withCMake-quickins) Note 2: **This is why Visual Studio 2019 is recommended above**, building HDF5 with newer versions of Visual Studio is not as clear (TODO: add VS2022 instructions). To clarify the HDF5Group setup/build instructions here is a short summary of the steps:
	1. Extract zip somewhere
	2. Go into `<path>\libraries\CMake-hdf5-1.10.8\CMake-hdf5-1.10.8` which contains .bat files etc and open command line (SHIFT + right click, open powershell window on windows)
	3. Execute the batch file or shell script for your platform. For example, `ctest -S HDF5config.cmake,BUILD_GENERATOR=VS201964 -C RelWithDebInfo -V -O hdf5.log` 
	4. Config `RelWithDebInfo` is the recommended configuration to allow running in Release mode or a debug mode (RelWithDebInfo).
	5. Wait for it to configure everything and run all its tests, most should pass.
	6. Locate and uncompress built binary its name is `HDF5-N-<platform>.<zip or tar.gz>` it appears in 2 locations, make sure to unzip the one in the build folder!
	7. Navagate to `<path>\libraries\CMake-hdf5-1.10.8\CMake-hdf5-1.10.8\build\HDF5-1.10.8-win64\HDF5-1.10.8-win64\share\cmake` or similar location to the cmake folder. Check that `hdf5-config.cmake` files are here.
	8. In general, to use HDF5 library users must first set the `HDF5_DIR` environment variable to the installed location of the CMake configuration files for HDF5. For example, on Windows the following path might be set:
		* `HDF5_DIR=C:\Users\USER\code\libraries\CMake-hdf5-1.10.8\CMake-hdf5-1.10.8\build\HDF5-1.10.8-win64\HDF5-1.10.8-win64\share\cmake`

3. Recommended (For use with Python): Install and build PyChrono according to [instructions.](https://api.projectchrono.org/module_python_installation.html) Note, instructions for building PyChrono include installing SWIG to generate python files-this is used again to generate python from HydroChrono C++ code.

* (optional) Gnuplot [Gnuplot Home](http://www.gnuplot.info/) or other plotting software

## Building HydroChrono (and demos)
1. Install and build above requirements, exact versions listed are required for proper functionality.
2. Now to build HydroChrono library: 
	1. Clone this project into a directory, and set up a build folder (i.e. clone project into HydroChrono directory and set up HydroChrono_build adjacent directory). 
	2. Just like when building Project Chrono, open CMake GUI specifying the location of the source files and built binaries for HydroChrono (HydroChrono and HydroChrono_build respectively). Configure and generate the solution setting 
		* `Chrono_DIR` as Chrono Build location (`../chrono_build/cmake`)
		* HDF5_DIR as (`../CMake-hdf5-1.10.8/CMake-hdf5-1.10.8/build/HDF5-1.10.8-win64/HDF5-1.10.8-win64/share/cmake` or similar path to find the `cmake` file at the end of this path). Note: version 1.10.8 of HDF5 works best with Visual Studio 2019.
		* Enable `HYDROCHRONO_ENABLE_DEMOS`, `HYDROCHRONO_ENABLE_IRRLICHT`, and `HYDROCHRONO_ENABLE_TESTS` to enable each feature. it is recommended to enable all of these, and the Irrlicht module depends on having Project Chrono be built with Irrlicht module enabled.
	3. Navigate to the build folder and open the solution in Visual Studio (or simply press "Open Project" in CMake GUI). Build the solution for HydroChrono in RelWithDebInfo mode (The `ALL_BUILD` project is the best for building and linking everything).
3. From Project Chrono build directory copy `chrono_build/bin/data` file into `HydroChrono_build/data` for optional shaders and logos
4. Navigate to `chrono_build/bin/RelWithDebInfo` folder and copy all .dll and .pdb files (not for demos) and paste them into `HydroChrono_build/demos/RelWithDebInfo` file. List of all files to copy:
	* ChronoEngine.dll
	* ChronoEngine.pdb
	* ChronoEngine_irrlicht.dll
	* ChronoEngine_irrlicht.pdb
	* ChronoEngine_postprocess.dll
	* ChronoEngine_postprocess.pdb
	* Irrlicht.dll (no pdb for Irrlicht)
	* and any other chrono module ChronoEngine_moduleName.dll and ChronoEngine_moduleName.pdb
5. Navigate to `HydroChrono_build/demos/RelWithDebInfo` and run demos.
	* When running demos, a command line argument is needed for the location of input files. To set this, run demo executables from the command line with the command argument for the absolute file location for `<path>/HydroChrono/demos` (Ex: the command `> ./sphere_decay.exe C:\Users\USERNAME\HydroChrono\demos` when executed from the `HydroChrono_build/demos/RelWithDebInfo` file will run the sphere heave decay test demo with the proper input file locations.)
		* Alternatively, set the command line argument for debug mode in a demo's properties. Do this in Visual Studio 2019 by navigating to the solution explorer, right clicking a demo, select `properties > Configuration Properties > Debugging` and set the option for Command Arguments to `C:\Users\USERNAME\HydroChrono\demos`.
	* Optionally copy plot files into results files and generate plots.

## Tests

If `HYDROCHRONO_ENABLE_TESTS` is enabled during configuration (need python) you can launch every tests in the build directory with the command:

	ctest .

or for launch only a subset of the tests (demos, examples, small, medium, long):

	ctest -L demos -L medium





## Files
* src - contains source code for HydroChrono library
	* hydro_forces.cpp
		* Defines classes to calculate hydrodynamic forces and apply to a multibody system.
	* chloadaddedmass.cpp
		* Defines a class to apply the added mass at infinite frequency to a system using Project Chrono Loadable objects.
	* h5fileinfo.cpp
		* Defines a class to read an h5 file and store system properties.
	* helper.cpp
		* Defines helper functions for various operations.
* include/hydroc - contains header files for above corresponding .cpp files.
* demos - contains information for various demos
	* demos/rm3 - contains required files for rm3 demos (decay test and regular wave test)
		* rm3_decay.cpp and rm3_reg_waves.cpp - demos showing rm3 decay and regular wave test respectively
		* demos/rm3/geometry - .obj files defining verticies on rm3 objects
		* demos/rm3/hydroData - contains rm3.h5 file
		* demos/rm3/postprocessing - contains .plt plotting files and WEC-Sim comparison data for plots. Copy these into results folder to plot demo outputs
	* demos/sphere - contains required files for sphere demos (decay test and regular wave test)
		* sphere_decay.cpp and sphere_reg_waves.cpp - demos showing sphere decay and regular wave test respectively
		* demos/sphere/geometry - .obj files defining verticies on sphere object
		* demos/sphere/hydroData - contains sphere.h5 file
		* demos/sphere/postprocessing - contains .plt and .py plotting files and WEC-Sim comparison data for plots. Copy these into results folder to plot demo outputs
	* demos/f3of - contains required files for f3of demo (decay test and regular wave test)
		* F3OF.cpp - demo showing f3of decay, optionally add regular waves
		* demos/f3of/geometry - .obj files defining verticies on f3of objects
		* demos/f3of/hydroData - contains f3of.h5 file
		* demos/f3of/postprocessing - contains .plt plotting files and WEC-Sim comparison data for plots. Copy these into results folder to plot demo outputs
	* demos/oswec - contains required files for oswec demo (decay test)
		* oswec_decay.cpp - demo showing oswec decay and optionally regular waves
		* demos/oswec/geometry - .obj files defining verticies on oswec object
		* demos/oswec/hydroData - contains oswec.h5 file
		* demos/oswec/postprocessing - contains .plt and .py plotting files and WEC-Sim comparison data for plots. Copy these into results folder to plot demo outputs
	* demos/Test - can be ignored
	* demos/python - can be ignored for now
* tests - contains multiple base tests used when option `HYDROCHRONO_ENABLE_TESTS` is enabled