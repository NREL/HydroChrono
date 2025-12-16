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

# Import the comparison template
sys.path.append(os.path.join(os.path.dirname(__file__), '../utilities'))
from compare_template import run_consistency_comparison

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
    safe_test_name = test_name.lower().replace(' ', '_').replace('-', '_')
    print(f"Plot filename: {plots_dir}/{safe_test_name}_comparison.png")
    y_label = "Heave (m)"
    tolerance = 1e-6
    
    # Run the consistency comparison using the template
    max_diff, passed = run_consistency_comparison(
        fname_rst, test_name, y_label,
        tolerance=tolerance,
        ref_label="Spectrum",
        sim_label="Eta-Import"
    )
    
    sys.exit(0 if passed else 1)
