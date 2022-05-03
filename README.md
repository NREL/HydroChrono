# Implementing Hydro Forces with Project Chrono and HDF5 files

## Requirements
* Chrono installed and built (currently using dev branch from Feb. 23, 2022) and Chrono dependencies (see [Chrono Install/Build Guide](https://api.projectchrono.org/tutorial_install_chrono.html) )
	* C++ Compiler
	* Eigen3
	* Irrlicht Engine
	* CMake
	* GIT client
* Install H5Cpp header file from HDF5Group [Download Link](https://portal.hdfgroup.org/display/support/Downloads) [Instructions](https://portal.hdfgroup.org/display/support/Building+HDF5+with+CMake#BuildingHDF5withCMake-quickins) Note: Version 1.10.8 or any 1.10.x recommended, other versions may not work as intended.
* Recommended (For use with Python): Install and build PyChrono according to [instructions.](https://api.projectchrono.org/module_python_installation.html) Note, instructions for building PyChrono include installing SWIG to generate python files-this is used again to generate python from HydroChrono C++ code.
* (optional) Gnuplot [Gnuplot Home](http://www.gnuplot.info/) or other plotting software

## Building HydroChrono (and demos)
1. Install and build above requirements.
2. Now to build HydroChrono library. 
	1. Clone this project into a directory, and set up a build folder (ie clone project into HydroChrono directory and set up HydroChrono_build adjacent directory). 
	2. Just like when building Project Chrono, open CMake GUI specifying the location of the source files and built binaries (HydroChrono and HydroChrono_build respectively). Configure and generate the solution, specifying Chrono_DIR as Chrono Build location (`../chrono_build/cmake`) and HDF5_DIR as (`../CMake-hdf5-1.10.8/CMake-hdf5-1.10.8/build/HDF5-1.10.8-win64/HDF5-1.10.8-win64/share/cmake`). Note: the current version of HDF5 works best with Visual Studio 2019 (hDF5Group offers alternatives for other versions)
	3. Navigate to the build folder and open the solution in Visual Studio. Build the solution for HydroChrono in Release mode (for debug mode see instructions for installing other configurations of HDF5 library) (The ALL_BUILD project is the best for building and linking everything).
3. From Project Chrono build directory copy `chrono_build/bin/data` file into `HydroChrono_build/data` for shaders and logos
4. Navigate to `HydroChrono_build/Release` (or Debug if set up properly) and run executables (Important note: may need to set up output files `../HydroChrono_build/Release/outfile/output.txt`)

## Files
* hydro_forces.cpp and hydro_forces.h
	* header and implementation files for hydro forces initialized through H5 files
* sphere_decay_demo.cpp
	* demo for hyrdo forces 
* sphere.h5 
	* file in project folder to define sphere data