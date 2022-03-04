#Implementing Hydro Forces with Project Chrono and HDF5 files

##Requirements
* Chrono Installed (currently using dev branch from Feb. 23, 2022) and Chrono dependencies (see [Chrono Install Guide](https://api.projectchrono.org/tutorial_install_chrono.html) )
	* C++ Compiler
	* Eigen3
	* Irrlicht Engine
	* CMake
	* GIT client
* H5Cpp header file from HDF5Group [github link](https://github.com/steven-varga/h5cpp)
* (optional) Gnuplot [Gnuplot Home](http://www.gnuplot.info/) or other plotting software

##Building Demos
1. Install above requirements, and build Chrono
2. Configure and generate CMake, specifying Chrono_DIR as Chrono Build location (`../chrono_build/cmake`) and HDF5_DIR as (`../CMake-hdf5-1.10.8/CMake-hdf5-1.10.8/build/HDF5-1.10.8-win64/HDF5-1.10.8-win64/share/cmake)
3. Build `test_for_chrono` in adjacent file `test_for_chrono_build`
4. Navigate to `test_for_chrono_build/Release` (or Debug) and run executables

##Files
* H5_force_classes.cpp and H5_force_classes.h
	* header and implementation files for hydro forces initialized through H5 files
* sphere_decay_test.cpp
	* demo for hyrdo forces 
* addedMassTest.cpp
	* demo for hydro forces with added mass implemented as custom load
	* addedMass class declared/defined in this file for now- eventually will be moved to H5_force_class