# Manual Diagnostic Tests

This directory contains diagnostic tools for **visual verification** of the state-space radiation damping implementation in HydroChrono. These are NOT automated tests - they generate comparison figures showing convolution vs state-space methods.

## Quick Start

```bash
# Build
cmake --build build --config Release --target sphere_decay_ss_comparison

# Run (generates CSV data files)
cd build/bin/tests/manual
./sphere_decay_ss_comparison

# Plot (generates verification figure)
python plot_model_comparison.py sphere --save
```

---

## Available Diagnostics

### Sphere Decay Test

Single-body heave decay simulation comparing convolution and state-space radiation.

```bash
cmake --build build --config Release --target sphere_decay_ss_comparison
./sphere_decay_ss_comparison
python plot_model_comparison.py sphere --save
```

**Output**: `output/sphere_ss_verification.png`

### OSWEC Decay Test

Two-body pitch decay simulation (flap + base) with increased fitting parameters for complex coupled kernels.

```bash
cmake --build build --config Release --target oswec_decay_ss_comparison
./oswec_decay_ss_comparison
python plot_model_comparison.py oswec --save
```

**Output**: `output/oswec_ss_verification.png`

---

## Output Files

### Data Files (CSV)

| File | Description |
|------|-------------|
| `{model}_kernel_fit.csv` | RIRF kernel vs state-space fit for key DOF |
| `{model}_convolution.csv` | Time series from convolution method |
| `{model}_state_space.csv` | Time series from state-space method |
| `{model}_performance.csv` | Computation time scaling data |

### Figures (PNG)

Each verification figure (`{model}_ss_verification.png`) contains:

1. **Top row**: Kernel fit quality (RIRF vs SS fit for key DOF)
2. **Middle row**: Time series comparison (convolution vs state-space)
3. **Bottom row**: Performance scaling (total time, per-step cost, speedup)

---

## Directory Structure

```
tests/manual/
├── sphere_decay_ss_comparison.cpp  # Sphere heave decay test
├── oswec_decay_ss_comparison.cpp   # OSWEC pitch decay test
├── plot_model_comparison.py        # Unified plotting script
├── CMakeLists.txt
└── README.md

build/bin/tests/manual/
├── sphere_decay_ss_comparison.exe
├── oswec_decay_ss_comparison.exe
├── plot_model_comparison.py
└── output/
    ├── sphere_ss_verification.png
    ├── oswec_ss_verification.png
    └── *.csv (data files)
```

---

## What These Diagnostics Show

### Kernel Fit Quality

Shows how accurately the Hankel-SVD fitter approximates the RIRF kernel:
- Blue line: Actual RIRF from BEM data
- Red dashed: State-space reconstruction
- Good fit: Lines overlap, R² close to 1.0

### Time Series Comparison

Verifies the state-space model produces forces equivalent to direct convolution:
- Blue: Convolution-based simulation
- Red dashed: State-space simulation
- Key metric: Maximum position/angle difference

### Performance Scaling

Demonstrates the computational advantage of state-space:
- **Convolution**: O(N) per step - time grows with simulation length
- **State-space**: O(1) per step - constant time regardless of history
- Speedup factor grows with simulation duration (typically 10-20× for long simulations)

---

## Requirements

Python dependencies:
```bash
pip install matplotlib pandas numpy
```

---

## Notes

- The sphere test uses default fitting parameters (`max_order=10`, `max_hankel_size=200`)
- The OSWEC test uses increased parameters (`max_order=15`, `max_hankel_size=500`) due to more complex coupled dynamics
- Both tests use the `RadiationStateSpaceComponent` integrated with `HydroSystem`
