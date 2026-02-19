#!/usr/bin/env python3
"""
OSWEC Irregular Waves Regression Test Comparison

Compares OSWEC irregular wave results (convolution radiation) against reference data.
"""

import sys
import os
from pathlib import Path

sys.path.append(os.path.join(os.path.dirname(__file__), '../utilities'))
from compare_template import run_comparison

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python compare_oswec_irreg_waves.py <reference_file> <results_file>")
        sys.exit(1)

    ref_file = sys.argv[1]
    results_file = sys.argv[2]
    print("Reference file: ", ref_file)
    print("Results file:   ", results_file)

    test_file_path = Path(results_file)
    plots_dir = test_file_path.parent / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    print(f"Plot will be saved to: {plots_dir}")

    test_name = "OSWEC Irregular Waves"
    y_label = "Flap Pitch (radians)"
    executable_patterns = ["oswec_irreg_waves_test", "oswec_irreg_waves"]
    pass_criteria = (1e-4, 0.02)

    n1, n2, passed = run_comparison(
        ref_file, results_file, test_name, y_label,
        executable_patterns, pass_criteria,
        status_name="oswec_irreg_waves"
    )

    sys.exit(0 if passed else 1)
