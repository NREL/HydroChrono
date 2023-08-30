.. _label-source_code_documentation:

Source Code Documentation
=========================

H5 Reader
---------

`h5fileinfo.cpp`

This file provides methods to manipulate and extract data from HDF5 files.

.. doxygenfile:: h5fileinfo.cpp
   :project: HydroChrono

Some Helper Functions
---------------------

`helper.cpp`

A collection of auxiliary functions that facilitate data manipulation and transformations.

.. doxygenfile:: helper.cpp
   :project: HydroChrono

Hydrodynamic Modeling
---------------------

`hydro_forces.cpp`
`hydro_forces.h`

This file encompasses the main hydrodynamic calculations and algorithms.

.. doxygenfile:: hydro_forces.cpp
   :project: HydroChrono

.. doxygenfile:: hydro_forces.h
   :project: HydroChrono


Wave Options
------------

`wave_types.cpp`

This source file includes various types of wave models, both regular and irregular.

.. doxygenfile:: wave_types.cpp
   :project: HydroChrono

Added Mass
----------

`chloadaddedmass.cpp`

Contains methods to include added mass in Chrono.

.. doxygenfile:: chloadaddedmass.cpp
   :project: HydroChrono