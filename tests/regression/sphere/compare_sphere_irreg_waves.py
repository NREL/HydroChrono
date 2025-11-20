#!/usr/bin/env python3
"""
Sphere Irregular Waves Regression Test Comparison Script

This script compares the results of the sphere irregular waves test against reference data.
"""

import sys
from pathlib import Path

# Add the utilities directory to the path to import the comparison template
sys.path.append(str(Path(__file__).parent.parent / "utilities"))
from compare_template import run_comparison


def main():
    """Main comparison function for sphere irregular waves test.

    When invoked from CTest, the reference and result file paths are passed as
    command-line arguments:
        compare_sphere_irreg_waves.py <ref_file> <result_file>

    If no arguments are provided (e.g., ad‑hoc runs), we fall back to a sane
    default based on the source and build tree layout.
    """

    # Prefer explicit paths from CTest if provided
    if len(sys.argv) == 3:
        ref_file = Path(sys.argv[1])
        result_file = Path(sys.argv[2])
    else:
        # Fallback: infer from repo layout (useful for manual runs)
        project_root = Path(__file__).resolve().parents[3]
        ref_file = project_root / "data" / "reference_data" / "sphere" / "hc_ref_sphere_irreg_waves.txt"
        result_file = (
            project_root / "build" / "bin" / "Release" / "results" / "tests" / "sphere" /
            "results_sphere_irreg_waves.txt"
        )

    if not result_file.exists():
        print(f"Error: Result file not found: {result_file}")
        sys.exit(1)

    print("Comparing sphere irregular waves test...")
    print(f"  Reference: {ref_file}")
    print(f"  Result:    {result_file}")

    # Run comparison using the template
    try:
        n1, n2, passed = run_comparison(
            str(ref_file),
            str(result_file),
            test_name="Sphere Irregular Waves",
            y_label="Heave (m)",
            executable_patterns=["sphere_irreg_waves_test"],
            pass_criteria=(1e-4, 0.02),
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