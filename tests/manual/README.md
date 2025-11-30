# Manual/Diagnostic Tests

This directory contains diagnostic tools for **visual verification** of HydroChrono components. These are NOT automated tests - they generate data files that you can plot and inspect.

## Available Diagnostics

### RadiationStateSpaceModel Verification

**Executable:** `radiation_ss_visual_check`  
**Python script:** `plot_radiation_ss.py`

Verifies the state-space radiation model by comparing numerical results against analytical solutions.

#### Scenarios covered:

1. **Constant velocity response** - Single mode, constant velocity input
   - Verifies exponential approach to steady state
   - Checks force = H × z relationship

2. **Step velocity response** - Velocity step change
   - Shows the "memory" effect of radiation damping
   - Tests transient response to input changes

3. **Multiple modes** - Superposition of slow and fast modes
   - Verifies that mode contributions add correctly
   - Shows different time constants interacting

4. **Time step sensitivity** - Same scenario with dt = 0.01, 0.1, 0.5, 1.0 s
   - Demonstrates unconditional stability of the exact integrator
   - Even very large time steps give reasonable results

## Usage

### 1. Build the diagnostic

```bash
cmake --build build --target radiation_ss_visual_check
```

### 2. Run it (generates CSV files)

```bash
# From the manual tests directory (recommended):
cd build/bin/tests/manual
./radiation_ss_visual_check

# This creates output/ subdirectory with CSV files:
#   output/radiation_ss_constant_v.csv
#   output/radiation_ss_step_v.csv
#   ...

# Or specify a custom output directory:
./radiation_ss_visual_check my_output_dir
```

### 3. Plot the results

```bash
# From the manual tests directory:
cd build/bin/tests/manual
python plot_radiation_ss.py

# This reads from output/ and displays plots

# Save plots to PNG files:
python plot_radiation_ss.py --save

# Or specify a custom data directory:
python plot_radiation_ss.py my_output_dir --save
```

## Directory Structure

After running, you'll have:
```
build/bin/tests/
├── unit/
│   ├── test_radiation_ss_model.exe
│   └── test_hydro_yaml_parser.exe
├── manual/
│   ├── radiation_ss_visual_check.exe
│   ├── plot_radiation_ss.py
│   └── output/
│       ├── radiation_ss_constant_v.csv
│       ├── radiation_ss_step_v.csv
│       ├── radiation_ss_multi_mode.csv
│       ├── radiation_ss_dt_sensitivity.csv
│       ├── plot_constant_velocity.png  (if --save)
│       ├── plot_step_velocity.png
│       ├── plot_multi_mode.png
│       └── plot_dt_sensitivity.png
└── regression/
    └── ...
```

## Requirements

The plotting script requires:
- Python 3.6+
- matplotlib
- pandas

Install with:
```bash
pip install matplotlib pandas
```

## Output Files

| File | Description |
|------|-------------|
| `radiation_ss_constant_v.csv` | Constant velocity response data |
| `radiation_ss_step_v.csv` | Step velocity response data |
| `radiation_ss_multi_mode.csv` | Multiple modes superposition data |
| `radiation_ss_dt_sensitivity.csv` | Time step sensitivity data |

## Why Visual Checks?

Unit tests verify *numerical correctness* - the numbers match expected values within tolerance. But visual checks help to build intuition about how the model behaves, spot unexpected patterns and communicate results to others effectively.
