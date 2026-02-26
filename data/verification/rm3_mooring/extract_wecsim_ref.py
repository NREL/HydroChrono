#!/usr/bin/env python3
"""
Extract reference data from a WEC-Sim/MoorDyn co-simulation of the RM3 mooring case.

Reference case: WEC-Sim Applications Mooring
https://github.com/WEC-Sim/WEC-Sim_Applications/tree/2d7569271c00e8736c2755192ee45941ffd49403/Mooring

Run that case to obtain the .mat outputs (e.g. etaData.mat and
Mooring/Mooring/MoorDyn/output/RM3MoorDyn_matlabWorkspace.mat). Then run this
script with --wecsim-dir pointing to the MoorDyn case directory (the directory
containing etaData.mat and the output/ subdirectory).

Reads the WEC-Sim MATLAB workspace (.mat) and wave elevation file, and writes
plain-text reference files suitable for HydroChrono verification tests.

Output files (committed to repo):
    data/verification/rm3_mooring/reference/wecsim_moordyn_body_motions.txt
    data/verification/rm3_mooring/reference/wecsim_moordyn_fairlead_tensions.txt
    data/verification/rm3_mooring/reference/wecsim_moordyn_wave_elevation.txt
    data/verification/rm3_mooring/reference/README.md
    data/verification/rm3_mooring/inputs/eta_rm3_mooring.txt

Usage:
    python data/verification/rm3_mooring/extract_wecsim_ref.py --wecsim-dir PATH

    --wecsim-dir is required: path to the WEC-Sim MoorDyn case directory.
"""

import argparse
import os
import sys
from pathlib import Path
from datetime import datetime

import numpy as np

try:
    import scipy.io as sio
except ImportError:
    sys.exit("scipy is required: pip install scipy")

try:
    import h5py
except ImportError:
    sys.exit("h5py is required: pip install h5py")


def extract(wecsim_dir: Path, repo_root: Path):
    eta_mat_path = wecsim_dir / "etaData.mat"
    workspace_mat_path = wecsim_dir / "output" / "RM3MoorDyn_matlabWorkspace.mat"

    for p in (eta_mat_path, workspace_mat_path):
        if not p.exists():
            sys.exit(f"Source file not found: {p}")

    ref_dir = repo_root / "data" / "verification" / "rm3_mooring" / "reference"
    inp_dir = repo_root / "data" / "verification" / "rm3_mooring" / "inputs"
    ref_dir.mkdir(parents=True, exist_ok=True)
    inp_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Wave elevation (etaData.mat -> HydroChrono eta format) ────────
    print(f"Reading {eta_mat_path} ...")
    eta_raw = sio.loadmat(str(eta_mat_path))

    eta_key = [k for k in eta_raw if not k.startswith("__")]
    if len(eta_key) != 1:
        sys.exit(f"Expected one variable in etaData.mat, found: {eta_key}")
    eta_data = eta_raw[eta_key[0]]  # (N, 2): [time, elevation]
    print(f"  eta shape: {eta_data.shape}, t=[{eta_data[0,0]:.2f}, {eta_data[-1,0]:.2f}]s")

    # HydroChrono eta format: "time : elevation" per line
    eta_out = inp_dir / "eta_rm3_mooring.txt"
    with open(eta_out, "w") as f:
        for row in eta_data:
            f.write(f"{row[0]:.6f} : {row[1]:.10e}\n")
    print(f"  -> {eta_out}  ({len(eta_data)} points)")

    # ── 2. WEC-Sim/MoorDyn workspace (.mat v7.3 via h5py) ───────────────
    print(f"Reading {workspace_mat_path} ...")
    ws = h5py.File(str(workspace_mat_path), "r")
    refs = ws["#refs#"]

    time = refs["Fb"]["time"][0, :]          # (40001,)
    elevation = refs["Fb"]["elevation"][0, :] # (40001,)
    n_steps = len(time)
    dt = time[1] - time[0]
    print(f"  time: {n_steps} steps, dt={dt:.4f}s, t=[{time[0]:.2f}, {time[-1]:.2f}]s")

    # Body positions: refs/8 = float (body1), refs/9 = plate (body2)
    # Shape (6, 40001): rows are [x, y, z, rx, ry, rz]
    float_pos = refs["8"][:]  # body1 (float)
    plate_pos = refs["9"][:]  # body2 (plate)
    print(f"  float Z at t=0: {float_pos[2, 0]:.4f}m  (expect -0.72)")
    print(f"  plate Z at t=0: {plate_pos[2, 0]:.4f}m  (expect -21.50)")

    # Fairlead tensions
    ft4 = refs["Cb"]["Lines"]["FairTen4"][0, :]
    ft5 = refs["Cb"]["Lines"]["FairTen5"][0, :]
    ft6 = refs["Cb"]["Lines"]["FairTen6"][0, :]
    print(f"  FairTen4 range: [{ft4.min():.1f}, {ft4.max():.1f}] N")

    ws.close()

    # ── 3. Write reference body motions ──────────────────────────────────
    body_out = ref_dir / "wecsim_moordyn_body_motions.txt"
    header = "Time(s)  FloatHeaveZ(m)  PlateHeaveZ(m)"
    body_arr = np.column_stack([time, float_pos[2, :], plate_pos[2, :]])
    np.savetxt(str(body_out), body_arr, header=header, fmt="%.6f", comments="")
    print(f"  -> {body_out}  ({n_steps} rows)")

    # ── 4. Write reference fairlead tensions ─────────────────────────────
    tension_out = ref_dir / "wecsim_moordyn_fairlead_tensions.txt"
    header = "Time(s)  FairTen4(N)  FairTen5(N)  FairTen6(N)"
    tension_arr = np.column_stack([time, ft4, ft5, ft6])
    np.savetxt(str(tension_out), tension_arr, header=header, fmt="%.6f", comments="")
    print(f"  -> {tension_out}  ({n_steps} rows)")

    # ── 5. Write reference wave elevation ────────────────────────────────
    wave_out = ref_dir / "wecsim_moordyn_wave_elevation.txt"
    header = "Time(s)  Elevation(m)"
    wave_arr = np.column_stack([time, elevation])
    np.savetxt(str(wave_out), wave_arr, header=header, fmt="%.6f", comments="")
    print(f"  -> {wave_out}  ({n_steps} rows)")

    # ── 6. Write README ──────────────────────────────────────────────────
    readme_out = ref_dir / "README.md"
    readme_out.write_text(f"""\
# WEC-Sim/MoorDyn RM3 Mooring Reference Data

## Provenance

Extracted on **{datetime.now().strftime('%Y-%m-%d %H:%M')}** by
`data/verification/rm3_mooring/extract_wecsim_ref.py` from a WEC-Sim/MoorDyn
co-simulation of the RM3 device with 6-line catenary mooring.

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
| MoorDyn          | 3 lines, 6 segments (anchor–clump + clump–fairlead each); input: data/yaml/rm3/mooring/lines_rm3.txt |

## Files

| File | Description |
|------|-------------|
| `wecsim_moordyn_body_motions.txt` | Time, float heave (Z), plate heave (Z) from WEC-Sim |
| `wecsim_moordyn_fairlead_tensions.txt` | Time, FairTen4, FairTen5, FairTen6 from MoorDyn (inside WEC-Sim) |
| `wecsim_moordyn_wave_elevation.txt` | Time, wave elevation as used in the co-simulation |
""", encoding="utf-8")
    print(f"  -> {readme_out}")

    print("\nDone.")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--wecsim-dir", type=str, required=True,
                        help="Path to WEC-Sim MoorDyn case directory (contains etaData.mat and output/)")
    args = parser.parse_args()

    wecsim_dir = Path(args.wecsim_dir)
    repo_root = Path(__file__).resolve().parent.parent

    extract(wecsim_dir, repo_root)


if __name__ == "__main__":
    main()
