#!/usr/bin/env python3
"""
OSWEC Regular Waves Regression Test Comparison Script

This script compares the results of the OSWEC regular waves test against reference data.
The test runs 16 different wave conditions (periods 4-24s), so we need to compare each one.
"""

import sys
import os
import glob
from pathlib import Path

# Add the utilities directory to the path to import the comparison template
sys.path.append(str(Path(__file__).parent.parent / "utilities"))
from compare_template import run_comparison, run_multi_column_comparison

def main():
    """Main comparison function for OSWEC regular waves test."""
    
    if len(sys.argv) != 3:
        print("Usage: python compare_oswec_reg_waves.py <reference_file> <test_file>")
        sys.exit(1)

    # Get reference and results files
    ref_file = sys.argv[1]
    results_file = sys.argv[2]
    print("Reference file: ", ref_file)
    print("Results file:   ", results_file)

    # Get the reference and results directories
    ht = os.path.split(ref_file)
    ref_dir = Path(ht[0])
    if not ref_dir.exists():
        print(f"Error: Reference directory not found: {ref_dir}")
        sys.exit(1)
    
    ht = os.path.split(results_file)
    results_dir = Path(ht[0])
    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}")
        sys.exit(1)
    
    # Find all result files produced by the C++ test
    # C++ writes files named: results_oswec_reg_waves_<N>.txt
    result_files = list(results_dir.glob("results_oswec_reg_waves_*.txt"))
    # Filter out duration files
    result_files = [f for f in result_files if "_duration" not in f.stem]
    result_files.sort(key=lambda x: int(x.stem.split('_')[-1]))
    
    if not result_files:
        print(f"Error: No result files found in {results_dir}")
        sys.exit(1)
    
    print(f"Found {len(result_files)} result files in {results_dir}")
    
    # Compare each wave condition
    all_passed = True
    
    for result_file in result_files:
        # Extract wave number from filename
        wave_num = result_file.stem.split('_')[-1]
        ref_file = ref_dir / f"hc_ref_oswec_reg_waves_{wave_num}.txt"
        
        if not ref_file.exists():
            print(f"Warning: Reference file {ref_file} not found, skipping wave {wave_num}")
            continue
        
        print(f"\nComparing wave condition {wave_num}...")
        print(f"  Reference: {ref_file}")
        print(f"  Result:    {result_file}")
        
        # Run comparison using the template
        try:
            n1, n2, passed = run_comparison(
                str(ref_file),
                str(result_file),
                test_name=f"OSWEC Regular Waves - Wave {wave_num}",
                y_label="Pitch (radians)",
                executable_patterns=["test_oswec_reg_waves"],
                pass_criteria=(1e-4, 0.02)
            )
            
            if not passed:
                all_passed = False
                print(f"  FAILED Wave {wave_num} comparison")
            else:
                print(f"  PASSED Wave {wave_num} comparison")
                
        except Exception as e:
            print(f"  ERROR comparing wave {wave_num}: {e}")
            all_passed = False
    
    if all_passed:
        print("\nAll wave conditions passed comparison!")
        sys.exit(0)
    else:
        print("\nSome wave conditions failed comparison!")
        sys.exit(1)

if __name__ == "__main__":
    main()
