#######################
User guide
#######################

.. todo::
	This is an example. Should be replaced

Global parameters
=================


:math:`Q_{\lambda}^{\rm ini} = \frac{F}{3600} \sum_{k=1}^{N_{\lambda}} V_{\lambda,k}\, c_{s,\lambda,k,\rm ini}`,

.. csv-table:: global parameters & variables
	:file: glb_parameters.csv
	:header-rows: 1
	:widths: 30, 12, 12, 12, 50

	

Input files
============

.. rubric:: Command file

The main input file (or *command* file) of HydroChrono is a json file, which can be named arbitrarily (e.g., 'cmd.json'
in the example below). It should contain a list of JSON objects with the following information:

	-  "Mesh": type and name of input file describing the lithium cell geometry and the simulation grid
	-  "Materials": name of input file providing all material/interface properties
	-  "InitialCondition": name of input file defining the initial condition
	-  "Runspec": name of input file for run specifications (cell operation, timestep control, output data,...)

.. code-block:: json
	:caption: JSON command file example
	:name: json-cmd-file

	{
		"Mesh":
		{
			"Type": "Structured-VTK",
			"Filename": "geometry.vts"
		},
		"Materials":
		{
			"Filename": "materials.json"
		},
		"InitialCondition":
		{
			"Filename": "initialcondition.json"
		},
		"Runspec":
		{
			"Filename": "runspec.json"
		}
	}


.. rubric:: Mesh file
