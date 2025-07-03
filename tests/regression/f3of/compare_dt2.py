#!/usr/bin/env python3
"""
F3OF DT2 Regression Test Comparison Script

This script compares HydroChrono results for F3OF DT2 (pitch decay test) 
against reference data from the demos folder.
"""

import sys
import os
from pathlib import Path
import numpy as np

# Import the comparison template
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from compare_template import create_comparison_plot

def main():
    # F3OF DT2 specific configuration
    test_name = "F3OF DT2 Pitch Decay"
    executable_patterns = ["f3of_dt2_test", "f3of_dt2_test.exe"]
    
    # Data files - use absolute paths
    hc_data_file = "results/CHRONO_F3OF_DT2_PITCH.txt"
    ref_data_file = os.path.join(os.path.dirname(__file__), "..", "reference_data", "f3of", "dt2", "hc_ref_f3of_dt2_pitch.txt")
    
    # Check if files exist
    if not os.path.exists(hc_data_file):
        print(f"Error: Test data file not found: {hc_data_file}")
        sys.exit(1)
    
    if not os.path.exists(ref_data_file):
        print(f"Error: Reference data file not found: {ref_data_file}")
        sys.exit(1)
    
    try:
        # Load data for validation
        ref_data = np.loadtxt(ref_data_file, skiprows=1)
        test_data = np.loadtxt(hc_data_file, skiprows=1)
        
        # Extract pitch column (column 2) for comparison
        ref_pitch_data = np.column_stack((ref_data[:, 0], ref_data[:, 2]))  # time, pitch
        test_pitch_data = np.column_stack((test_data[:, 0], test_data[:, 2]))  # time, pitch
        
        # Show where the plots will be saved
        test_file_path = Path(hc_data_file)
        plots_dir = test_file_path.parent / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        print(f"Plots will be saved to: {plots_dir}")
        
        # Find the executable path
        executable_path = None
        if executable_patterns:
            from compare_template import find_executable
            executable_path = find_executable(test_file_path.parent, executable_patterns)
        
        # Generate comparison plot
        def rel_to_root(path):
            try:
                project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
                return os.path.relpath(path, project_root)
            except Exception:
                return str(path)
        
        n1, n2 = create_comparison_plot(
            ref_pitch_data, test_pitch_data, test_name, plots_dir, 
            ref_file_path=rel_to_root(ref_data_file), 
            test_file_path=rel_to_root(hc_data_file),
            executable_path=rel_to_root(str(executable_path)) if executable_path else None,
            y_label="Base Pitch (rad)",
            executable_patterns=executable_patterns
        )
        
        # Additional F3OF-specific validations
        
        # Check base surge (column 1) - should be zero for DT2
        diff_base_surge = np.linalg.norm(ref_data[:,1] - test_data[:,1]) / len(ref_data[:,1])
        if diff_base_surge > 1e-10:
            print(f"F3OF DT2 validation failed: Base surge difference {diff_base_surge:.2e} > 1e-10")
            sys.exit(1)
        
        # Check base pitch (column 2) - main test for DT2
        diff_base_pitch = np.linalg.norm(ref_data[:,2] - test_data[:,2]) / len(ref_data[:,2])
        if diff_base_pitch > 1e-6:
            print(f"F3OF DT2 validation failed: Base pitch difference {diff_base_pitch:.2e} > 1e-6")
            sys.exit(1)
        
        # Check flap fore pitch (column 3) - should match base pitch for DT2
        diff_flap_fore_pitch = np.linalg.norm(ref_data[:,3] - test_data[:,3]) / len(ref_data[:,3])
        if diff_flap_fore_pitch > 1e-6:
            print(f"F3OF DT2 validation failed: Flap fore pitch difference {diff_flap_fore_pitch:.2e} > 1e-6")
            sys.exit(1)
        
        # Check flap aft pitch (column 4) - should match base pitch for DT2
        diff_flap_aft_pitch = np.linalg.norm(ref_data[:,4] - test_data[:,4]) / len(ref_data[:,4])
        if diff_flap_aft_pitch > 1e-6:
            print(f"F3OF DT2 validation failed: Flap aft pitch difference {diff_flap_aft_pitch:.2e} > 1e-6")
            sys.exit(1)
        
        # Check template comparison results
        l2_threshold, linf_threshold = 1e-6, 1e-6
        if (n1 > l2_threshold or n2 > linf_threshold):
            print(f"F3OF DT2 TEST FAILED - L2 Norm: {n1:.2e}, L-infinity Norm: {n2:.2e}")
            sys.exit(1)
        else:
            print(f"F3OF DT2 TEST PASSED - L2 Norm: {n1:.2e}, L-infinity Norm: {n2:.2e}")
            print(f"Generated plot: {plots_dir}/{test_name.replace(' ', '_')}_comparison.png")
        
    except Exception as e:
        print(f"Error during comparison: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 