#!/usr/bin/env python3
"""
Sphere Irregular Waves Regression Test Comparison Script

This script compares the results of the sphere irregular waves test against reference data.
"""

import sys
import os
from pathlib import Path

# Add the utilities directory to the path to import the comparison template
sys.path.append(str(Path(__file__).parent.parent / "utilities"))
from compare_template import run_comparison


def main():
    """Main comparison function for sphere irregular waves test."""

    if len(sys.argv) != 3:
        print("Usage: python compare.py <reference_file> <test_file>")
        sys.exit(1)

    # Get reference and results files
    ref_file = sys.argv[1]
    results_file = sys.argv[2]
    print("Reference file: ", ref_file)
    print("Results file:   ", results_file)

    if not os.path.exists(results_file):
        print(f"Error: Result file not found: {results_file}")
        sys.exit(1)

    print("Comparing sphere irregular waves test...")
    print(f"  Reference: {ref_file}")
    print(f"  Result:    {results_file}")

    # Run comparison using the template
    try:
        n1, n2, passed = run_comparison(
            ref_file,
            results_file,
            test_name="Sphere Irregular Waves",
            y_label="Heave (m)",
            executable_patterns=["sphere_irreg_waves_test"],
            pass_criteria=(1e-4, 0.02),
            status_name="sphere_irreg_waves"
        )

        if passed:
            print("Sphere irregular waves comparison PASSED")
            sys.exit(0)
        else:
            print("Sphere irregular waves comparison FAILED")
            sys.exit(1)

    except Exception as e:
        print(f"ERROR during comparison: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()