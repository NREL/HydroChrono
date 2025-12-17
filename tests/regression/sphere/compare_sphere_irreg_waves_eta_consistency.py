#!/usr/bin/env python3
"""
HydroChrono Sphere Irregular Waves Eta Consistency Test Comparison

This script validates that spectrum-generated and eta-imported irregular waves
produce identical results using the standardized comparison template.

Usage:
    python compare_sphere_irreg_waves_eta_consistency.py <reference_file> <results_file>
    Note: reference_file is unused for this self-validating test
"""

import sys
import os
from pathlib import Path
import numpy as np

# Import the comparison template
sys.path.append(os.path.join(os.path.dirname(__file__), '../utilities'))
from compare_template import create_comparison_plot, format_path

if __name__ == '__main__':
    """
    Compare spectrum-generated vs eta-imported irregular waves for consistency
    """
    if len(sys.argv) != 3:
        print("Usage: python compare_sphere_irreg_waves_eta_consistency.py <reference_file> <results_file>")
        print("Note: reference_file is unused for this self-validating test")
        sys.exit(1)

    # Reference file is unused - this is a self-validating test
    fname_rst = sys.argv[2]

    if not os.path.exists(fname_rst):
        print(f"ERROR: Results file not found: {fname_rst}")
        sys.exit(1)

    # Show where the plot will be saved
    test_file_path = Path(fname_rst)
    plots_dir = test_file_path.parent / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    print(f"Plot will be saved to: {plots_dir}")

    # Test-specific configuration
    test_name = "Sphere Irreg Waves Eta Consistency"
    y_label = "Heave (m)"
    tolerance = 1e-6

    # Load the results file (columns: time, heave_spectrum, heave_eta, diff)
    # Skip header lines (title + separator)
    data = np.loadtxt(fname_rst, skiprows=2)
    
    # Extract columns as ref (spectrum) and test (eta) data for comparison
    ref_data = np.column_stack((data[:, 0], data[:, 1]))   # time, spectrum
    test_data = np.column_stack((data[:, 0], data[:, 2]))  # time, eta
    
    print(f"Loaded {len(data)} data points")

    # Generate comparison plot using standard template
    n1, n2 = create_comparison_plot(
        ref_data, test_data, test_name, plots_dir,
        ref_file_path=format_path(fname_rst),
        test_file_path=format_path(fname_rst),
        y_label=y_label
    )

    # Check pass/fail
    max_diff = np.max(np.abs(data[:, 3]))  # Use pre-computed diff column
    passed = max_diff < tolerance
    
    if passed:
        print(f"TEST PASSED - Max difference: {max_diff:.2e}, Tolerance: {tolerance:.2e}")
    else:
        print(f"TEST FAILED - Max difference: {max_diff:.2e} exceeds tolerance: {tolerance:.2e}")
    
    sys.exit(0 if passed else 1)
