---
layout: page
title: Sphere Model (IEA OES Task 10) - Verification
parent_section: Verification
---

## Overview

Developed under the IEA OES Task 10 project, the Sphere Model serves as a simplified wave energy converter (WEC) for verifying numerical modeling software. This verification encompasses various tests, including still water decay and dynamic responses in regular and irregular wave conditions, some incorporating a damping coefficient to represent power take-off (PTO).

<p align="center">
  <img src="{{ site.baseurl }}/verification/images_sphere/sphere_mesh_and_irreg.png" alt="Visualization of the Sphere. Left: the mesh used to compute hydrodynamic coefficients. Right: visualization of the sphere simulated with HydroChrono in irregular waves." width="75%" />
</p>

## Model Parameters

Key physical parameters of the Sphere Model, like center of buoyancy, center of gravity, and volume displacement, are essential for hydrodynamic analysis:

| Name                        | Symbol                         | Value                  | Units       |
| --------------------------- | ------------------------------ | ---------------------- | ----------- |
| Water Density               | \\( \rho_{\text{water}} \\)    | \\( 1 \times 10^3 \\)  | kg/m³       |
| Gravity                     | \\( g \\)                      | 9.81                   | m/s²        |
| Water Depth                 | \\( d \\)                      | \\( \infty \\)         | m           |
| Sphere Mass                 | \\( m \\)                      | \\( 261.8 \times 10^3 \\) | kg       |
| Sphere Radius               | \\( r \\)                      | 5                      | m           |
| Sphere Centre               | \\( \mathbf{C}_{\text{centre}} \\) | \\( [0.0, 0.0, 0.0] \\) | m         |
| Sphere Center of Gravity    | \\( \mathbf{C}_{\text{G}} \\)  | \\( [0.0, 0.0, -2.0] \\) | m          |

## Results

The model's heave response data is analyzed for verification:

<p align="center">
  <img src="{{ site.baseurl }}/verification/images_sphere/sphere_decay_1m_verification.png" alt="Sphere decay test verification results (comparisons against selected IEA OES Task 10 participants)." width="75%" />
</p>

The verification results establish the model's reliability in simulating hydrodynamic behaviors. While these tests offer a strong foundation for confidence in the model's performance, ongoing verification efforts are necessary to continually verify and validate the code.

## References

For detailed information about the Sphere Model, refer to the following source:

- Wendt, F. F., et al. (2017). International Energy Agency Ocean Energy Systems Task 10 Wave Energy Converter Modeling Verification and Validation. Retrieved from [OSTI.gov](https://www.osti.gov/biblio/1401957)
