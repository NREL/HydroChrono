#!/usr/bin/env python3
"""
F3OF DT3 Regression Test Comparison Script

This script compares HydroChrono results for F3OF DT3 (flap pitch decay test)
against reference data from the demos folder.
"""

import sys
import os
from pathlib import Path
import numpy as np

# Import the comparison template
sys.path.append(os.path.join(os.path.dirname(__file__), '../utilities'))
from compare_template import run_multi_column_comparison, clip_to_common_time

def main():
    # F3OF DT3 specific configuration
    test_name = "F3OF DT3 Flap Pitch Decay"
    executable_patterns = ["f3of_dt3_test", "f3of_dt3_test.exe"]

    if len(sys.argv) != 3:
        print("Usage: python compare.py <reference_file> <test_file>")
        sys.exit(1)

    # Get reference and results files
    ref_file = sys.argv[1]
    results_file = sys.argv[2]
    print("Reference file: ", ref_file)
    print("Results file:   ", results_file)

    # Check if files exist
    if not os.path.exists(results_file):
        print(f"Error: Test data file not found: {results_file}")
        sys.exit(1)

    if not os.path.exists(ref_file):
        print(f"Error: Reference data file not found: {ref_file}")
        sys.exit(1)

    try:
        # Load data for validation
        ref_data = np.loadtxt(ref_file, skiprows=1)
        test_data = np.loadtxt(results_file, skiprows=1)
        ref_data, test_data = clip_to_common_time(ref_data, test_data)

        # Show where the plots will be saved
        test_file_path = Path(results_file)
        plots_dir = test_file_path.parent / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        print(f"Plots will be saved to: {plots_dir}")

        # Find the executable path
        executable_path = None
        if executable_patterns:
            from compare_template import find_executable
            executable_path = find_executable(test_file_path.parent, executable_patterns)

        # Define the columns to plot and their configurations
        test_configs = [
            {
                'column_index': 3,  # Flap Fore Pitch
                'test_name': f"{test_name} - Flap Fore Pitch",
                'y_label': "Flap Fore Pitch (radians)",
                'validation_tolerance': (1e-6, 1e-6),
                'status_name': 'f3of_dt3_flap_fore_pitch'
            },
            {
                'column_index': 4,  # Flap Aft Pitch
                'test_name': f"{test_name} - Flap Aft Pitch",
                'y_label': "Flap Aft Pitch (radians)",
                'validation_tolerance': (1e-6, 1e-6),
                'status_name': 'f3of_dt3_flap_aft_pitch'
            }
        ]

        # Run the multi-column comparison using the template
        results = run_multi_column_comparison(
            ref_file, results_file, test_configs,
            executable_patterns=executable_patterns
        )

        # Additional F3OF-specific validations (from the original script)
        # Interpolate reference onto test time grid to handle different durations
        nval = test_data.shape[0]
        x = np.linspace(test_data[0, 0], test_data[nval-1, 0], nval)

        # Check base surge (column 1)
        ref_col1 = np.interp(x, ref_data[:,0], ref_data[:,1])
        test_col1 = np.interp(x, test_data[:,0], test_data[:,1])
        diff_base_surge = np.linalg.norm(ref_col1 - test_col1) / nval
        if diff_base_surge > 1e-6:
            print(f"F3OF DT3 validation failed: Base surge difference {diff_base_surge:.2e} > 1e-6")
            sys.exit(1)

        # Check base pitch (column 2) - more stringent tolerance
        ref_col2 = np.interp(x, ref_data[:,0], ref_data[:,2])
        test_col2 = np.interp(x, test_data[:,0], test_data[:,2])
        diff_base_pitch = np.linalg.norm(ref_col2 - test_col2) / nval
        if diff_base_pitch > 1e-10:
            print(f"F3OF DT3 validation failed: Base pitch difference {diff_base_pitch:.2e} > 1e-10")
            sys.exit(1)

        # Check flap fore pitch (column 3)
        ref_col3 = np.interp(x, ref_data[:,0], ref_data[:,3])
        test_col3 = np.interp(x, test_data[:,0], test_data[:,3])
        diff_flap_fore_pitch = np.linalg.norm(ref_col3 - test_col3) / nval
        if diff_flap_fore_pitch > 1e-6:
            print(f"F3OF DT3 validation failed: Flap fore pitch difference {diff_flap_fore_pitch:.2e} > 1e-6")
            sys.exit(1)

        # Check flap aft pitch (column 4)
        ref_col4 = np.interp(x, ref_data[:,0], ref_data[:,4])
        test_col4 = np.interp(x, test_data[:,0], test_data[:,4])
        diff_flap_aft_pitch = np.linalg.norm(ref_col4 - test_col4) / nval
        if diff_flap_aft_pitch > 1e-6:
            print(f"F3OF DT3 validation failed: Flap aft pitch difference {diff_flap_aft_pitch:.2e} > 1e-6")
            sys.exit(1)

        # Check if all template comparisons passed
        all_passed = all(result[2] for result in results)

        if all_passed:
            print("F3OF DT3 TEST PASSED - All comparisons within tolerance")
            print(f"Generated plots:")
            for config in test_configs:
                plot_name = config['test_name'].lower().replace(' ', '_').replace('-', '_')
                print(f"  - {plots_dir}/{plot_name}_comparison.png")
        else:
            print("F3OF DT3 TEST FAILED - Some comparisons outside tolerance")
            sys.exit(1)

    except Exception as e:
        print(f"Error during comparison: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 