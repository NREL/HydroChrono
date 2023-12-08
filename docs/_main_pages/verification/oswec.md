---
layout: page
title: Oscillating Surge Wave Energy Converter (OSWEC) - Verification
parent_section: Verification
---

## Overview

The Oscillating Surge Wave Energy Converter (OSWEC) is the verification case tested in this section, using the example model provided in the WEC-Sim package ([WEC-Sim Reference](https://wec-sim.github.io/WEC-Sim/dev/user/tutorials.html#oscillating-surge-wec-oswec)). The verification process serves as a comparative measure against WEC-Sim results for assessing the accuracy of HydroChrono.

<p align="center">
  <img src="{{ site.baseurl }}/verification/images_oswec/oswec_model.png" alt="Visualization of OSWEC. Left: hydrodynamic mesh for coefficient computation. Right: OSWEC system in the Chrono GUI." width="50%" />
</p>




## Model Parameters

The OSWEC model's fundamental physical parameters—such as hinge location, center of gravity, and flap mass—are shown in the following table.

| Name                         | Symbol                        | Value           | Units         |
| ---------------------------- | ----------------------------- | --------------- | ------------- |
| Hinge Location, \\( x_{Hinge} \\) | \\( [0.0, 0.0, -8.9] \\)      | -               | m             |
| Flap Center of Gravity, \\( C_{G} \\) | \\( [0.0, 0.0, -3.9] \\)  | -               | m             |
| Flap Mass, \\( m_{flap} \\)      | \\( 127 \times 10^3 \\)      | -               | kg            |
| Flap Pitch Inertia, \\( I_{yy, flap} \\) | \\( 1.85 \times 10^6 \\) | -               | kg.m^2        |
| Water Depth, \\( d_{water} \\)    | \\( 10.9 \\)                 | -               | m             |
| Water Density, \\( \rho_{water} \\) | \\( 1000 \\)               | -               | kg/m^3        |


## Results

The OSWEC pitch decay tests and Response Amplitude Operators (RAOs) were focused on for HydroChrono verification.

<p align="center">
  <img src="{{ site.baseurl }}/verification/images_oswec/oswec_pitch_decay_verification.png" alt="OSWEC 10-degree pitch decay test, displaying a high degree of correlation between HydroChrono and WEC-Sim results." width="50%" />
</p>


<p align="center">
  <img src="{{ site.baseurl }}/verification/images_oswec/oswec_pitch_rao_verification.png" alt="RAOs from WEC-Sim and HydroChrono for OSWEC in regular waves, revealing minor discrepancies at the resonance frequency." width="50%" />
</p>


The verification of the OSWEC model provides confidence in the numerical implementation of the hydrodynamic forces for rotational degrees of freedom. The OSWEC's base is treated as a separate body in this system, providing some confidence in the multibody hydrodynamics implementation. However, since the base is fixed, further verification is required in order to check the accuracy of hydrodynamic interactions in multibody systems.

## References

For further information about the OSWEC model, please refer to the following:

- [WEC-Sim Tutorials: Oscillating Surge WEC (OSWEC)](https://wec-sim.github.io/WEC-Sim/dev/user/tutorials.html#oscillating-surge-wec-oswec)
- [WEC-Sim OSWEC Example](https://github.com/WEC-Sim/WEC-Sim/tree/master/examples/OSWEC)
