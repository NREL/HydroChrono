# WEC-Sim/MoorDyn RM3 Mooring Reference Data

## Provenance

Extracted on **2026-02-24 14:22** by
`data/verification/rm3_mooring/extract_wecsim_ref.py` from a WEC-Sim/MoorDyn
co-simulation of the RM3 device with a 3-line catenary mooring (6 MoorDyn segments: anchor–clump + clump–fairlead per line).

### Source

Reference data were produced by running the [WEC-Sim Applications Mooring](https://github.com/WEC-Sim/WEC-Sim_Applications/tree/2d7569271c00e8736c2755192ee45941ffd49403/Mooring) case (WEC-Sim/MoorDyn co-simulation of the RM3 device). That run generates MATLAB outputs (e.g. `etaData.mat` and `Mooring/Mooring/MoorDyn/output/RM3MoorDyn_matlabWorkspace.mat`). The plain-text reference files in this directory were then generated using `data/verification/rm3_mooring/extract_wecsim_ref.py` and are used by the HydroChrono verification tests and comparison plots.

### WEC-Sim simulation parameters

| Parameter        | Value                                           |
|------------------|-------------------------------------------------|
| dt               | 0.01 s                                          |
| End time         | 400 s                                           |
| Ramp time        | 40 s                                            |
| Wave type        | Elevation import (etaData.mat)                  |
| Water depth      | 70 m                                            |
| Float mass       | 725 834 kg (equilibrium)                        |
| Plate mass       | 886 691 kg (equilibrium)                        |
| Float inertia    | [20907301, 21306090.66, 37085481.11] kg-m^2     |
| Plate inertia    | [94419614.57, 94407091.24, 28542224.82] kg-m^2  |
| Plate init disp  | [0, 0, -0.21] m                                 |
| PTO damping     | 1 200 000 N-s/m                                 |
| PTO stiffness    | 0                                               |
| MoorDyn          | 3 lines, 6 segments (anchor–clump + clump–fairlead each); input: `data/yaml/rm3/mooring/lines_rm3.txt` |

## Files

| File | Description |
|------|-------------|
| `wecsim_moordyn_body_motions.txt` | Time, float heave (Z), plate heave (Z) from WEC-Sim |
| `wecsim_moordyn_fairlead_tensions.txt` | Time, FairTen4, FairTen5, FairTen6 from MoorDyn (inside WEC-Sim) |
| `wecsim_moordyn_wave_elevation.txt` | Time, wave elevation as used in the co-simulation |
