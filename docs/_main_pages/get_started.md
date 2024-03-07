---
layout: page
title: Get Started
---

## Getting Started with HydroChrono

Welcome to HydroChrono! To get started with using HydroChrono, follow the steps below to download and run the application.

### Downloading HydroChrono

1. Navigate to the [Releases](https://github.com/NREL/HydroChrono/releases) page of the HydroChrono repository.
2. Look for the latest release (currently v0.3) and download the `hydrochrono_wsi.exe` file.

### Running HydroChrono

To run HydroChrono, you'll need to use the command line. Here's a basic overview of the usage:

`hydrochrono_wsi.exe [InputFile.m] [resultsDir] [--gui/--nogui] [--quiet] [--nohydro]`


- `[InputFile.m]`: The path to the input file (in WEC-Sim format) that contains your simulation configurations.
- `[resultsDir]`: The directory where the simulation results will be saved.
- `--gui/--nogui`: Include `--gui` to run the simulation with a graphical user interface or `--nogui` to run it without.
- `--quiet`: Include this option if you want to suppress detailed logging output.
- `--nohydro`: Include this option if you want to disable hydrodynamic calculations.

### Example Usage

In the `HydroChrono\examples\wsi` directory, you'll find example .m files that you can use to test the application. Here's how you might run a simulation with one of these example files:

`hydrochrono_wsi.exe HydroChrono\examples\wsi\rm3\wecSimInputFile.m results --gui`

This command runs a simulation using `rm3\wecSimInputFile.m` as the input file, saves the results to a directory called `results`, and runs without the GUI.

### Warning

Please be aware that HydroChrono is still under development, and not all WEC-Sim functionalities are available through the .m input file. Users may encounter bugs or unexpected behaviors.

### Troubleshooting & Reporting Issues

If you encounter any issues while using HydroChrono's examples, please check the following:

- Ensure you've entered the command-line arguments correctly.
- Verify that the input file path and results directory exist and are accessible.
- Check the console output for any error messages that might indicate what went wrong.

If you encounter any issues or have suggestions for improvement, please don't hesitate to report them on the [Issues page](https://github.com/NREL/HydroChrono/issues) of our GitHub repository. Your feedback is invaluable to us and helps make HydroChrono better for everyone!

Thank you for using HydroChrono!
