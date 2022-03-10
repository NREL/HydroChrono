# Implementing Hydro Forces with Project Chrono and HDF5 files

## Requirements
* Chrono Installed (currently using dev branch from Feb. 23, 2022) and Chrono dependencies (see [Chrono Install Guide](https://api.projectchrono.org/tutorial_install_chrono.html) )
	* C++ Compiler
	* Eigen3
	* Irrlicht Engine
	* CMake
	* GIT client
* H5Cpp header file from HDF5Group [github link](https://github.com/steven-varga/h5cpp)
* (optional) Gnuplot [Gnuplot Home](http://www.gnuplot.info/) or other plotting software

## Building Demo
1. Install above requirements, and build Chrono
2. Configure and generate CMake, specifying Chrono_DIR as Chrono Build location (`../chrono_build/cmake`) and HDF5_DIR as (`../CMake-hdf5-1.10.8/CMake-hdf5-1.10.8/build/HDF5-1.10.8-win64/HDF5-1.10.8-win64/share/cmake`)
3. Build `test_for_chrono` in adjacent file `test_for_chrono_build`
4. Copy `chrono_build/bin/data` file into `test_for_chrono_build/data`
5. Navigate to `test_for_chrono_build/Release` (or Debug) and run executables (Important note: may need to set up output files `../test_for_chrono_build/Release/outfile/output.txt`)

## Files
* hydro_forces.cpp and hydro_forces.h
	* header and implementation files for hydro forces initialized through H5 files
* sphere_decay_demo.cpp
	* demo for hyrdo forces 
* sphere.h5 
	* file in project folder to define sphere data