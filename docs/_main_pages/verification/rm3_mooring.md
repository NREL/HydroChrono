---
layout: page
title: RM3 Mooring (WEC-Sim/MoorDyn) - Verification
parent_section: Verification
---

## Overview

This verification case compares HydroChrono with MoorDyn against the WEC-Sim/MoorDyn co-simulation of the RM3 (Reference Model 3) two-body point absorber with a 3-line catenary mooring (each line has two segments—anchor to clump weight, then clump to fairlead—giving 6 MoorDyn line segments in total). The reference case is the [WEC-Sim Applications Mooring](https://github.com/WEC-Sim/WEC-Sim_Applications/tree/2d7569271c00e8736c2755192ee45941ffd49403/Mooring) example. Running that case produces MATLAB outputs that were converted to plain-text reference data using `data/verification/rm3_mooring/extract_wecsim_ref.py`; those reference files live in `data/verification/rm3_mooring/reference/` and are used by the HydroChrono verification tests and comparison plots.

<p align="center">
  <img src="{{ site.baseurl }}/verification/images_rm3_mooring/rm3_mooring_screenshot.png" alt="RM3 floating platform with mooring lines colored by tension (HydroChrono/MoorDyn simulation)." width="75%" />
</p>

## Model and input parameters

Key parameters match the WEC-Sim/MoorDyn RM3 mooring setup. The verification test uses the same timestep, ramp, wave elevation import, body masses and inertias, PTO settings, and MoorDyn configuration as the reference.

| Parameter        | Value                                           |
|------------------|-------------------------------------------------|
| Timestep         | 0.01 s                                          |
| Simulation       | 60 s comparison window (after 40 s ramp)        |
| Ramp time        | 40 s                                            |
| Wave type        | Elevation import (irregular, from reference)    |
| Water depth      | 70 m                                            |
| Float mass       | 725 834 kg                                      |
| Plate mass       | 886 691 kg                                      |
| Float inertia    | [20907301, 21306090.66, 37085481.11] kg·m²      |
| Plate inertia    | [94419614.57, 94407091.24, 28542224.82] kg·m²   |
| Plate init disp  | [0, 0, -0.21] m                                 |
| PTO damping      | 1 200 000 N·s/m                                 |
| PTO stiffness    | 0                                               |
| MoorDyn          | 3 lines, 6 segments (anchor–clump + clump–fairlead each); see below |

Reference data were produced from the [WEC-Sim Applications Mooring](https://github.com/WEC-Sim/WEC-Sim_Applications/tree/2d7569271c00e8736c2755192ee45941ffd49403/Mooring) case and extracted with `data/verification/rm3_mooring/extract_wecsim_ref.py`.

### Mooring configuration

The mooring is defined in the MoorDyn input file `data/yaml/rm3/mooring/lines_rm3.txt`. There are 3 physical mooring lines, each composed of two chain segments: anchor → clump weight (240 m) and clump weight → fairlead on the plate (40 m). MoorDyn treats these as 6 line segments; the verification compares fairlead tensions for the three lines attached to the body (FairTen4, FairTen5, FairTen6). The MoorDyn input file contents are:

```
---------------------- LINE TYPES -----------------------------------------------------
LineType  Diam      Mass/m        EA          BA/-zeta   EI         Cd     Ca    CdAx     CaAx    
(-)       (m)       (kg/m)        (N)          (N-s/-)  (N-m^2)    (-)     (-)    (-)     (-) 
Chain      0.144      126.0     583.376E6       -0.8      0        1.6     1.0    0.05     0.0
---------------------------- BODIES -----------------------------------------------------
ID   Attachment  X0     Y0    Z0     r0      p0     y0     Mass  CG*   I*      Volume   CdA*   Ca
(#)     (-)      (m)    (m)   (m)   (deg)   (deg)  (deg)   (kg)  (m)  (kg-m^2)  (m^3)   (m^2)  (-)
1       Coupled     0    0     -21.5    0       0      0     0.0   0    0         0       0      0
---------------------- POINTS  -----------------------------------------------------
ID      Attachment    X      Y     Z            Mass   Volume     CdA     CA
(#)     (word/ID)    (m)    (m)   (m)           (kg)   (mˆ3)      (m^2)   (-)
1         Body1     -3.0     0      11.50       0        0         0      0
2         Fixed     -267.0   0      -70.00       0        0         0      0
3         Body1     1.5     2.598   11.50       0        0         0      0
4         Fixed     133.5   231.23  -70.00       0        0         0      0
5         Body1     1.5    -2.598   11.50       0        0         0      0
6         Fixed     133.5  -231.23  -70.00       0        0         0      0
7         Free      -40.0     0     -10.00     16755   33.510     12.566   1
8         Free      20.0   34.642   -10.00     16755   33.510     12.566   1
9         Free      20.0  -34.642   -10.00     16755   33.510     12.566   1
---------------------- LINES -----------------------------------------------------
ID     LineType  AttachA  AttachB  UnstrLen     NumSegs   LineOutputs
(#)      (-)       (#)      (#)       (m)         (-)       (-)
1        Chain     2         7        240.0        15        tp
2        Chain     4         8        240.0        15        tp
3        Chain     6         9        240.0        15        tp
4        Chain     7         1         40.0         5        tp
5        Chain     8         3         40.0         5        tp
6        Chain     9         5         40.0         5        tp
---------------------- SOLVER OPTIONS-----------------------------------------
0.0005   dtM          - time step to use in mooring integration
0        WaveKin      - wave kinematics flag (0=neglect, the only option currently supported)
3.0e6    kBot         - bottom stiffness
3.0e5    cBot         - bottom damping
70       WtrDpth      - water depth
5.0      CdScaleIC    - factor by which to scale drag coefficients during dynamic relaxation IC gen
0.001    threshIC     - threshold for IC con
200.0    TmaxIC       - max time for ic gen (s)
2        writeLog     - Write a log file
0        WriteUnits   - 0: do not write the units header on the output files
1        disableOutTime - 1: do not write timesteps to command window
-------------------------- OUTPUTS --------------------------------
FairTen4
FairTen5
FairTen6
--------------------- need this line ------------------
```

Lines 1–3 run from anchor points (2, 4, 6) to clump weights (7, 8, 9); lines 4–6 run from clump weights to the plate (body) attachment points (1, 3, 5). FairTen4, FairTen5, FairTen6 are the tensions at the fairlead end of lines 4, 5, 6.

## Verification results

The verification compares body heave (float and plate) and fairlead tensions (FairTen4, FairTen5, FairTen6) from HydroChrono/MoorDyn against the WEC-Sim/MoorDyn reference. Below is the FairTen4 (fourth fairlead tension) comparison: the top panel shows reference vs simulation time series; the bottom panel shows the error (reference minus simulation).

<p align="center">
  <img src="{{ site.baseurl }}/verification/images_rm3_mooring/rm3_mooring_verification_fairten4_comparison.png" alt="RM3 Mooring verification: FairTen4 reference vs simulation and error (Ref − Sim)." width="85%" />
</p>

Typical metrics for this run: correlation ≈ 0.998, mean error ≈ 46 N, max error ≈ 3 kN, L∞ norm ≈ 3 kN. Mean fairlead tension is on the order of 170 kN, so the max error is about 1.8% of mean tension.

## Discrepancies

Remaining differences are small and are likely due to solver and integrator choices: HydroChrono uses Chrono’s HHT (Hilber–Hughes–Taylor) implicit integrator with GMRES, while WEC-Sim may use different time integration (e.g. Runge–Kutta). Small differences can also come from MoorDyn coupling details and wave-elevation interpolation. The close agreement (high correlation, low relative error) supports confidence in the coupled HydroChrono/MoorDyn implementation for this configuration.

## References

- [WEC-Sim Applications – Mooring](https://github.com/WEC-Sim/WEC-Sim_Applications/tree/2d7569271c00e8736c2755192ee45941ffd49403/Mooring) — reference case used to generate the RM3 mooring reference data
- Reference data and provenance: `data/verification/rm3_mooring/reference/README.md` in the HydroChrono repository
- MoorDyn input (mooring layout): `data/yaml/rm3/mooring/lines_rm3.txt`
