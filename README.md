# Implementing Hydro Forces with Project Chrono and HDF5 files

## Requirements
* Chrono installed and built in Release mode (currently using dev branch see note below) (built with Irrlicht module) and Chrono dependencies (see [Chrono Install/Build Guide](https://api.projectchrono.org/tutorial_install_chrono.html) )
	* C++ Compiler (Visual Studio 2019 or newer required, 2019 recommended)
	* Eigen3 (3.3.0 or newer required, 3.4.0 recommended)
	* Irrlicht Engine (1.8.4)
	* CMake
	* GIT client
* Note: to revert to the same commit of Chrono open Powershell in your Project Chrono directory and type `git reset --hard 0af74f59d`. The message `HEAD is now at 0af74f59d ...` confirms you are on the correct commit of Project Chrono dev branch.
* Install H5Cpp header file from HDF5Group [Download Link](https://portal.hdfgroup.org/display/support/Downloads) [Instructions](https://portal.hdfgroup.org/display/support/Building+HDF5+with+CMake#BuildingHDF5withCMake-quickins) Note: Version 1.10.8 or any 1.10.x recommended, other versions may not work as intended. Note 2: This is why Visual Studio 2019 is recommended above, building HDF5 with newer versions of Visual Studio is not as clear
* Recommended (For use with Python): Install and build PyChrono according to [instructions.](https://api.projectchrono.org/module_python_installation.html) Note, instructions for building PyChrono include installing SWIG to generate python files-this is used again to generate python from HydroChrono C++ code.
* (optional) Gnuplot [Gnuplot Home](http://www.gnuplot.info/) or other plotting software

## Building HydroChrono (and demos)
1. Install and build above requirements.
2. Now to build HydroChrono library. 
	1. Clone this project into a directory, and set up a build folder (ie clone project into HydroChrono directory and set up HydroChrono_build adjacent directory). 
	2. Just like when building Project Chrono, open CMake GUI specifying the location of the source files and built binaries (HydroChrono and HydroChrono_build respectively). Configure and generate the solution, specifying Chrono_DIR as Chrono Build location (`../chrono_build/cmake`) and HDF5_DIR as (`../CMake-hdf5-1.10.8/CMake-hdf5-1.10.8/build/HDF5-1.10.8-win64/HDF5-1.10.8-win64/share/cmake`). Note: version 1.10.8 of HDF5 works best with Visual Studio 2019 (hDF5Group offers alternatives for other versions)
	3. Navigate to the build folder and open the solution in Visual Studio. Build the solution for HydroChrono in Release mode (for debug mode see instructions for installing other configurations of HDF5 library) (The ALL_BUILD project is the best for building and linking everything).
3. From Project Chrono build directory copy `chrono_build/bin/data` file into `HydroChrono_build/data` for shaders and logos
4. Navigate to `HydroChrono_build/Release` (or Debug if set up properly) and run executables (Important note: may need to set up output files `../HydroChrono_build/Release/outfile/output.txt` and/or change file name for `sphere.h5`)

## Files
* hydro_forces.cpp and hydro_forces.h
	* header and implementation files for hydro forces initialized through H5 files
* sphere_decay_demo.cpp
	* demo for hyrdo forces 
* sphere.h5 
	* file in project folder to define sphere data