#!/usr/bin/env python3
"""
RM3 Regular Waves Regression Test Comparison

This script compares the RM3 regular waves test results with reference data and generates
comparison plots using the standardized template.
"""

import sys
import os
from pathlib import Path
import numpy as np

# Import the comparison template
sys.path.append(os.path.join(os.path.dirname(__file__), '../utilities'))
from compare_template import run_multi_column_comparison, write_status_file, clip_to_common_time

def main():

    if len(sys.argv) != 3:
        print("Usage: python compare.py <reference_file> <test_file>")
        sys.exit(1)

    # Get reference and results files
    ref_file = sys.argv[1]
    results_file = sys.argv[2]
    print("Reference file: ", ref_file)
    print("Results file:   ", results_file)
    
    # RM3 Regular Waves specific configuration
    test_name = "RM3 Regular Waves Test"
    executable_patterns = ["rm3_reg_waves_test", "rm3_reg_waves_test.exe"]
    
    # Define the columns to plot and their configurations
    test_configs = [
        {
            'column_index': 1,  # Float Heave
            'test_name': f"{test_name} - Float Heave",
            'y_label': "Float Heave (m)",
            'validation_tolerance': (1e-4, 0.02)  # RM3-specific tolerance
        },
        {
            'column_index': 2,  # Plate Heave
            'test_name': f"{test_name} - Plate Heave",
            'y_label': "Plate Heave (m)",
            'validation_tolerance': (1e-4, 0.02)  # RM3-specific tolerance
        },
        {
            'column_index': 3,  # Float Drift
            'test_name': f"{test_name} - Float Drift",
            'y_label': "Float Drift (m)",
            'validation_tolerance': (1e-4, 0.02)  # RM3-specific tolerance
        }
    ]
    
    try:
        # Load data for additional RM3-specific validations
        ref_data = np.loadtxt(ref_file, skiprows=1)  # Skip header line
        test_data = np.loadtxt(results_file, skiprows=1)  # Skip header line
        ref_data, test_data = clip_to_common_time(ref_data, test_data)
        
        # Show where the plots will be saved
        test_file_path = Path(results_file)
        plots_dir = test_file_path.parent / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        print(f"Plots will be saved to: {plots_dir}")

        # Create temporary files with interpolated data to match time steps
        import tempfile
        
        # Interpolate reference data to match test data time steps
        nval = test_data.shape[0]
        x = np.linspace(test_data[0, 0], test_data[nval-1, 0], nval)
        
        # Interpolate the float heave, plate heave, and float drift data
        floatHeaveRef = np.interp(x, ref_data[:,0], ref_data[:,1])
        plateHeaveRef = np.interp(x, ref_data[:,0], ref_data[:,2])
        floatDriftRef = np.interp(x, ref_data[:,0], ref_data[:,3])
        
        # Create interpolated reference data
        ref_data_interp = np.column_stack((x, floatHeaveRef, plateHeaveRef, floatDriftRef))
        
        # Create temporary files with interpolated data
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_ref:
            np.savetxt(temp_ref.name, ref_data_interp, fmt='%.6f')
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_test:
            np.savetxt(temp_test.name, test_data, fmt='%.6f')
        
        try:
            # Override the plots directory to ensure it's saved in the correct location
            import sys as sys_module
            sys_module.path.append(os.path.join(os.path.dirname(__file__), '..'))
            from compare_template import create_comparison_plot, format_path
            
            # Manually create the comparison plot to control the output directory
            for config in test_configs:
                column_index = config['column_index']
                test_name = config['test_name']
                y_label = config['y_label']
                
                # Create data arrays for this column
                ref_col_data = np.column_stack((ref_data_interp[:, 0], ref_data_interp[:, column_index]))
                test_col_data = np.column_stack((test_data[:, 0], test_data[:, column_index]))
                
                # Create the plot in the correct directory
                create_comparison_plot(
                    ref_col_data, test_col_data, test_name, plots_dir,
                    ref_file_path=format_path(ref_file),
                    test_file_path=format_path(results_file),
                    y_label=y_label,
                    executable_patterns=executable_patterns
                )
            
            # Create dummy results for compatibility
            results = [(0.0, 0.0, True)] * len(test_configs)
        finally:
            # Clean up temporary files
            os.unlink(temp_ref.name)
            os.unlink(temp_test.name)
        
        # Additional RM3-specific validations
        nval = test_data.shape[0]
        
        # Resample refData to testData sampling rate
        x = np.linspace(test_data[0, 0], test_data[nval-1, 0], nval)
        
        # Compare float heave
        floatHeaveRef = np.interp(x, ref_data[:,0], ref_data[:,1])
        floatHeaveTest = np.interp(x, test_data[:,0], test_data[:,1])
        floatHeaveComp = floatHeaveRef - floatHeaveTest

        # Compare plate heave
        plateHeaveRef = np.interp(x, ref_data[:,0], ref_data[:,2])
        plateHeaveTest = np.interp(x, test_data[:,0], test_data[:,2])
        plateHeaveComp = plateHeaveRef - plateHeaveTest

        # Compare float drift
        floatDriftRef = np.interp(x, ref_data[:,0], ref_data[:,3])
        floatDriftTest = np.interp(x, test_data[:,0], test_data[:,3])
        floatDriftComp = floatDriftRef - floatDriftTest

        # Frobenius norm - Float Heave
        floatHeaven1 = np.linalg.norm(floatHeaveComp)/nval
        # infinity norm - Float Heave
        floatHeaven2 = np.linalg.norm(floatHeaveComp, np.inf)

        # Frobenius norm - Plate Heave
        plateHeaven1 = np.linalg.norm(plateHeaveComp)/nval
        # infinity norm - Plate Heave
        plateHeaven2 = np.linalg.norm(plateHeaveComp, np.inf)

        # Frobenius norm - Float Drift
        floatDriftn1 = np.linalg.norm(floatDriftComp)/nval
        # infinity norm - Float Drift
        floatDriftn2 = np.linalg.norm(floatDriftComp, np.inf)
        
        if (floatHeaven1 > 1e-4 or floatHeaven2 > 0.02 or 
            plateHeaven1 > 1e-4 or plateHeaven2 > 0.02 or
            floatDriftn1 > 1e-4 or floatDriftn2 > 0.02):
            print(f"RM3 validation failed: Float heave difference {floatHeaven1:.2e} > 1e-4 or {floatHeaven2:.2e} > 0.02")
            print(f"RM3 validation failed: Plate heave difference {plateHeaven1:.2e} > 1e-4 or {plateHeaven2:.2e} > 0.02")
            print(f"RM3 validation failed: Float drift difference {floatDriftn1:.2e} > 1e-4 or {floatDriftn2:.2e} > 0.02")
            write_status_file(test_file_path.parent, "rm3_reg_waves", "FAIL")
            sys.exit(1)
        
        # Check if all template comparisons passed
        all_passed = all(result[2] for result in results)
        
        if all_passed:
            print("RM3 REGULAR WAVES TEST PASSED - All comparisons within tolerance")
            write_status_file(test_file_path.parent, "rm3_reg_waves", "PASS")
            print(f"Generated plots:")
            for config in test_configs:
                plot_name = config['test_name'].lower().replace(' ', '_').replace('-', '_')
                print(f"  - {config['test_name']}: {plots_dir}/{plot_name}_comparison.png")
        else:
            print("RM3 REGULAR WAVES TEST FAILED - Some comparisons outside tolerance")
            write_status_file(test_file_path.parent, "rm3_reg_waves", "FAIL")
            sys.exit(1)
        
    except Exception as e:
        print(f"Error during comparison: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 