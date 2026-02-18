#!/usr/bin/env python3
"""
HydroChrono Sphere Decay (State-Space Radiation) Regression Test Comparison

Compares sphere decay results using state-space radiation approximation against
the convolution-based reference data. Tolerances are slightly relaxed relative
to the convolution test since state-space is an O(1) approximation.
"""

import sys
import os
from pathlib import Path

sys.path.append(os.path.join(os.path.dirname(__file__), '../utilities'))
from compare_template import run_comparison

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python compare_sphere_decay_ss.py <reference_file> <results_file>")
        sys.exit(1)

    ref_file = sys.argv[1]
    results_file = sys.argv[2]
    print("Reference file: ", ref_file)
    print("Results file:   ", results_file)

    test_file_path = Path(results_file)
    plots_dir = test_file_path.parent / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    print(f"Plot will be saved to: {plots_dir}")

    test_name = "Sphere Decay Test (State-Space)"
    safe_test_name = test_name.lower().replace(' ', '_').replace('-', '_')
    print(f"Plot filename: {plots_dir}/{safe_test_name}_comparison.png")
    y_label = "Heave (m)"
    executable_patterns = ["sphere_decay_ss_test", "sphere_decay_ss"]
    pass_criteria = (5e-4, 0.05)

    n1, n2, passed = run_comparison(
        ref_file, results_file, test_name, y_label,
        executable_patterns, pass_criteria,
        status_name="sphere_decay_ss"
    )

    sys.exit(0 if passed else 1)
