# HydroChrono

**HydroChrono** is an emerging hydrodynamics simulation tool designed to model complex ocean systems. Seamlessly integrated with the [Project Chrono](https://projectchrono.org/) physics engine, it offers a powerful C++ API for a wide range of simulations.

## Capabilities

HydroChrono extends the versatility of Project Chrono to hydrodynamic applications, enabling you to simulate multibody wave energy converters (WECs) and floating offshore wind turbines (FOWTs) with Chrono::FEA for flexible bodies.

## Get Started

Visit the [HydroChrono website](https://nrel.github.io/HydroChrono/) for information about the underlying theory, build instructions, and usage details. Currently, HydroChrono is available as a source build, offering customization and optimization opportunities. A binary distribution with Python wrapper is in development and coming soon.

## Limitations

- Currently only supports first-order linear potential flow forces.
- Hydrodynamic forces are currently limited to rigid bodies.

## Vision and Future Goals

- Expand hydrodynamics to incorporate nonlinear and 2nd order forces.
- Integrate advanced simulation features using Project Chrono's DEM and FEA modules.
- Support seamless transitioning from potential flow to CFD and SPH for detailed FSI analysis.
- Develop a Python API to broaden accessibility and ease of use.
- Foster and support open-source collaboration in the research community.
