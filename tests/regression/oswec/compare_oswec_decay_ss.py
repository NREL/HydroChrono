#!/usr/bin/env python3
"""
OSWEC Decay (State-Space Radiation) Regression Test Comparison

Compares OSWEC decay results using state-space radiation approximation against
the convolution-based reference data. Tolerances are slightly relaxed relative
to the convolution test since state-space is an O(1) approximation.
"""

import sys
import os
from pathlib import Path
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '../utilities'))
from compare_template import run_multi_column_comparison, write_status_file, clip_to_common_time

def main():
    if len(sys.argv) != 3:
        print("Usage: python compare_oswec_decay_ss.py <reference_file> <results_file>")
        sys.exit(1)

    ref_file = sys.argv[1]
    results_file = sys.argv[2]
    print("Reference file: ", ref_file)
    print("Results file:   ", results_file)

    test_name = "OSWEC Decay Test (State-Space)"
    executable_patterns = ["oswec_decay_ss_test", "oswec_decay_ss_test.exe"]

    test_configs = [
        {
            'column_index': 1,
            'test_name': f"{test_name} - Flap Pitch",
            'y_label': "Flap Pitch (radians)",
            'validation_tolerance': (5e-4, 0.05)
        }
    ]

    try:
        ref_data = np.loadtxt(ref_file, skiprows=1)
        test_data = np.loadtxt(results_file, skiprows=1)
        ref_data, test_data = clip_to_common_time(ref_data, test_data)

        test_file_path = Path(results_file)
        plots_dir = test_file_path.parent / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        print(f"Plots will be saved to: {plots_dir}")

        import tempfile

        nval = test_data.shape[0]
        x = np.linspace(test_data[0, 0], test_data[nval-1, 0], nval)

        flapPitchRef = np.interp(x, ref_data[:,0], ref_data[:,1])
        ref_data_interp = np.column_stack((x, flapPitchRef))

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_ref:
            np.savetxt(temp_ref.name, ref_data_interp, fmt='%.6f')

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_test:
            test_data_interp = np.column_stack((test_data[:,0], test_data[:,1]))
            np.savetxt(temp_test.name, test_data_interp, fmt='%.6f')

        try:
            import sys as sys_module
            sys_module.path.append(os.path.join(os.path.dirname(__file__), '..'))
            from compare_template import create_comparison_plot, format_path

            for config in test_configs:
                column_index = config['column_index']
                config_test_name = config['test_name']
                y_label = config['y_label']

                ref_col_data = np.column_stack((ref_data_interp[:, 0], ref_data_interp[:, column_index]))
                test_col_data = np.column_stack((test_data_interp[:, 0], test_data_interp[:, column_index]))

                create_comparison_plot(
                    ref_col_data, test_col_data, config_test_name, plots_dir,
                    ref_file_path=format_path(ref_file),
                    test_file_path=format_path(results_file),
                    y_label=y_label,
                    executable_patterns=executable_patterns
                )

            results = [(0.0, 0.0, True)] * len(test_configs)
        finally:
            os.unlink(temp_ref.name)
            os.unlink(temp_test.name)

        nval = test_data.shape[0]
        x = np.linspace(test_data[0, 0], test_data[nval-1, 0], nval)
        flapPitchRef = np.interp(x, ref_data[:,0], ref_data[:,1])
        flapPitchTest = np.interp(x, test_data[:,0], test_data[:,1])
        flapPitchComp = flapPitchRef - flapPitchTest

        flapPitchn1 = np.linalg.norm(flapPitchComp)/nval
        flapPitchn2 = np.linalg.norm(flapPitchComp, np.inf)

        # Relaxed tolerances for state-space approximation
        if (flapPitchn1 > 5e-4 or flapPitchn2 > 0.05):
            print(f"OSWEC SS validation failed: Flap pitch difference {flapPitchn1:.2e} > 5e-4 or {flapPitchn2:.2e} > 0.05")
            write_status_file(test_file_path.parent, "oswec_decay_ss", "FAIL")
            sys.exit(1)

        all_passed = all(result[2] for result in results)

        if all_passed:
            print("OSWEC DECAY (STATE-SPACE) TEST PASSED - All comparisons within tolerance")
            write_status_file(test_file_path.parent, "oswec_decay_ss", "PASS")
            print(f"Generated plots:")
            for config in test_configs:
                plot_name = config['test_name'].lower().replace(' ', '_').replace('-', '_')
                print(f"  - {config['test_name']}: {plots_dir}/{plot_name}_comparison.png")
        else:
            print("OSWEC DECAY (STATE-SPACE) TEST FAILED - Some comparisons outside tolerance")
            write_status_file(test_file_path.parent, "oswec_decay_ss", "FAIL")
            sys.exit(1)

    except Exception as e:
        print(f"Error during comparison: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
