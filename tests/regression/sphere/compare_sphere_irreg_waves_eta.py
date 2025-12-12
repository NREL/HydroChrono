#!/usr/bin/env python3
"""
Sphere Irregular Waves with ETA Import Regression Test Comparison Script

This script compares the results of the sphere irregular waves with ETA import test against reference data.
Note: This test uses the same reference data as the regular irregular waves test since the output format is identical.
"""

import sys
from pathlib import Path
import os

# Add the utilities directory to the path to import the comparison template
sys.path.append(str(Path(__file__).parent.parent.parent / "utilities"))
from compare_template import run_comparison

def main():
    """Main comparison function for sphere irregular waves with ETA import test."""
    
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
    
    # Run comparison using the template
    try:
        n1, n2, passed = run_comparison(
            ref_file,
            results_file,
            test_name="Sphere Irregular Waves with ETA Import",
            y_label="Heave (m)",
            executable_patterns=["sphere_irreg_waves_eta_test"],
            pass_criteria=(1e-4, 0.02)
        )
        
        if passed:
            print("Sphere irregular waves with ETA import comparison PASSED")
            sys.exit(0)
        else:
            print("Sphere irregular waves with ETA import comparison FAILED")
            sys.exit(1)
            
    except Exception as e:
        print(f"ERROR during comparison: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 