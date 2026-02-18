#!/usr/bin/env python3
"""
HydroChrono Sphere Irregular Waves with ETA Import Regression Test Comparison

This script compares sphere irregular waves with eta import test results against reference data
using the standardized comparison template.

Usage:
    python compare_sphere_irreg_waves_eta.py <reference_file> <test_file>
    python compare_sphere_irreg_waves_eta.py default <test_file>  # Uses default reference data
"""

import sys
import os
from pathlib import Path

# Import the comparison template
sys.path.append(os.path.join(os.path.dirname(__file__), '../utilities'))
from compare_template import run_comparison

if __name__ == '__main__':
    """
    Compare sphere irregular waves eta test results with reference data
    """

    if len(sys.argv) != 3:
        print("Usage: python compare.py <reference_file> <test_file>")
    sys.exit(1)

    # Get reference and results files
    ref_file = sys.argv[1]
    results_file = sys.argv[2]
    print("Reference file: ", ref_file)
    print("Results file:   ", results_file)

    # Show where the plot will be saved
    test_file_path = Path(results_file)
    plots_dir = test_file_path.parent / "plots"
    # Ensure the plots directory exists
    plots_dir.mkdir(parents=True, exist_ok=True)
    print(f"Plot will be saved to: {plots_dir}")

    # Test-specific configuration
    test_name = "Sphere Irregular Waves with ETA Import"
    # Convert test_name to lowercase with underscores for filename
    safe_test_name = test_name.lower().replace(' ', '_').replace('-', '_')
    print(f"Plot filename: {plots_dir}/{safe_test_name}_comparison.png")
    y_label = "Heave (m)"
    executable_patterns = ["test_sphere_irreg_waves_eta", "sphere_irreg_waves_eta"]
    pass_criteria = (1e-4, 0.02)  # (L2 threshold, L-infinity threshold)
    
    # Run the comparison using the template
    n1, n2, passed = run_comparison(
        ref_file, results_file, test_name, y_label, 
        executable_patterns, pass_criteria,
        status_name="sphere_irreg_waves_eta"
    )
    
    sys.exit(0 if passed else 1)
