.. _label-prerequisites:

Prerequisites
=============

Introduction
------------

This page provides a detailed list of prerequisites and installation guidelines for users who wish to use HydrChrono. It covers the necessary software, libraries, and additional tools.

Section 1: Chrono and Dependencies
----------------------------------

Before proceeding with the scientific software installation, users need to ensure they have Chrono v8.0.0 installed and built.

For detailed instructions on how to install and build Chrono, please refer to the `Chrono Install/Build Guide <https://api.projectchrono.org/tutorial_install_chrono.html>`_. For visualization, insure this build includes the Irrlicht module.

Below are the necessary dependencies and the recommended versions:

* C++ Compiler (Visual Studio 2019 or newer required, 2019 recommended)
* Eigen3 (3.3.0 or newer required, 3.4.0 recommended)
* Irrlicht Engine (1.8.4)
* CMake (3.18.2 recommended)
* GIT client


Section 2: H5Cpp from HDF5Group
-------------------------------

Users need to install the H5Cpp header file from HDF5Group. We recommend using version 1.10.8 as other versions may not be fully compatible. 

.. note::
   
   This is why Visual Studio 2019 is recommended. Building HDF5 with newer versions of Visual Studio may not be straightforward. Further instructions for building with VS2022 will be provided in the future.

To facilitate the installation of H5Cpp, here's a brief summary:

1. Download the 1.10.8 release .zip from `HDF5Group <https://portal.hdfgroup.org/display/support/Downloads>`_.
2. Extract the zip file to a preferred location.
3. Navigate to `<path>\\libraries\\CMake-hdf5-1.10.8\\CMake-hdf5-1.10.8` and open a command line or powershell window.
4. Run the respective batch file or shell script for your platform. For instance: `ctest -S HDF5config.cmake,BUILD_GENERATOR=VS201964 -C RelWithDebInfo -V -O hdf5.log`.
   
   - The `RelWithDebInfo` option is optional (users could choose `Debug` or `Release` etc.). However, it's important to maintain consistency with the configuration option used to build Chrono (be it `Release`, `Debug`, or `RelWithDebInfo`) and the one you plan to use for HydroChrono.
5. Follow through with the configuration process.
6. Once completed, locate the compressed binary named `HDF5-N-<platform>.<zip or tar.gz>` within the build folder and extract it.
7. Check the presence of `hdf5-config.cmake` files in the directory: `<path>\\libraries\\CMake-hdf5-1.10.8\\CMake-hdf5-1.10.8\\build\\HDF5-1.10.8-win64\\HDF5-1.10.8-win64\\share\\cmake`.
8. Users must set the `HDF5_DIR` environment variable to the installed location of the CMake configuration files for HDF5. For instance, on Windows: `HDF5_DIR=C:\\Users\\USER\\code\\libraries\\CMake-hdf5-1.10.8\\CMake-hdf5-1.10.8\\build\\HDF5-1.10.8-win64\\HDF5-1.10.8-win64\\share\\cmake`.

Section 4: Optional Tools
-------------------------

* Gnuplot or any other plotting software of your choice. Learn more at `Gnuplot Home <http://www.gnuplot.info/>`.
* For users interested in Python integration, it is recommended to install and build PyChrono. Follow the `instructions <https://api.projectchrono.org/module_python_installation.html>`_ provided by Project Chrono. This guide also covers the installation of SWIG for generating python files from C++ code.
