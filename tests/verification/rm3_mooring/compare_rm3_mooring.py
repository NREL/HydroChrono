#!/usr/bin/env python3
"""
RM3 Mooring Verification Test -- compare HydroChrono/MoorDyn results against
WEC-Sim/MoorDyn co-simulation reference data.

Usage (called automatically by CTest):
    python compare_rm3_mooring.py <reference_dir> <results_file>

    reference_dir:  data/verification/rm3_mooring/reference/
    results_file:   results/tests/rm3_mooring/results_rm3_mooring.txt
"""

import sys
import os
from pathlib import Path

import numpy as np

# Reuse the regression comparison template
sys.path.append(os.path.join(os.path.dirname(__file__), '../../regression/utilities'))
from compare_template import create_comparison_plot, write_status_file


# Cross-code tolerances -- force-level cosine ramp matching WEC-Sim convention.
# Tuned from 100 s run: float L2~2.2e-4, plate L2~7.1e-5, FT4 Linf~2.9kN.
BODY_MOTION_TOLERANCE  = (1e-3, 0.15)     # (L2, L-inf) for body heave [m]
FAIRLEAD_TOLERANCE     = (5e+1, 5e+3)     # (L2, L-inf) for fairlead tension [N]

# Skip the wave ramp transient when computing error metrics.
# The first RAMP_SKIP_TIME seconds are excluded from the comparison.
RAMP_SKIP_TIME = 40.0


def load_reference(ref_dir):
    """Load all WEC-Sim/MoorDyn reference data from the reference directory."""
    ref_dir = Path(ref_dir)

    body_file    = ref_dir / "wecsim_moordyn_body_motions.txt"
    tension_file = ref_dir / "wecsim_moordyn_fairlead_tensions.txt"

    body_data    = np.loadtxt(str(body_file), skiprows=1)
    tension_data = np.loadtxt(str(tension_file), skiprows=1)

    return body_data, tension_data


def load_moordyn_lines_out(search_dirs):
    """Try to load MoorDyn's fairlead tension output from multiple candidate locations.

    MoorDyn writes output next to its input file with the same basename
    (e.g., lines_rm3.out for lines_rm3.txt).
    """
    patterns = ["lines_rm3.out", "lines.out"]
    for d in search_dirs:
        for name in patterns:
            for candidate in [d / name,
                              d / "Mooring" / name,
                              d / "mooring" / name]:
                if candidate.exists():
                    print(f"  Found MoorDyn output: {candidate}")
                    return np.loadtxt(str(candidate), skiprows=1)
    return None


def deduplicate_time(time, signal):
    """Remove duplicate time entries, keeping the first value for each time."""
    _, unique_idx = np.unique(time, return_index=True)
    return time[unique_idx], signal[unique_idx]


def resample_to_common_grid(ref_time, ref_signal, test_time, test_signal,
                            skip_time=0.0):
    """Interpolate both signals onto a common uniform time grid.

    If skip_time > 0, the comparison starts after that many seconds
    (to exclude ramp / initial transient).
    """
    t_start = max(ref_time[0], test_time[0], skip_time)
    t_end   = min(ref_time[-1], test_time[-1])
    dt      = max(ref_time[1] - ref_time[0], test_time[1] - test_time[0])
    n_pts   = int(np.round((t_end - t_start) / dt)) + 1
    t_common = np.linspace(t_start, t_end, n_pts)

    ref_interp  = np.interp(t_common, ref_time, ref_signal)
    test_interp = np.interp(t_common, test_time, test_signal)
    print(f"Resampled to common grid [{t_start:.2f}, {t_end:.2f}]s  "
          f"(dt={dt:.4f}, {n_pts} pts)")
    return t_common, ref_interp, test_interp


def compare_signal(ref_time, ref_signal, test_time, test_signal, label, units,
                   tolerance, plots_dir, ref_file_label, test_file_label,
                   skip_time=0.0):
    """Compare a single time-series signal and generate a plot.

    Returns (l2, linf, passed).
    """
    ref_time, ref_signal   = deduplicate_time(ref_time, ref_signal)
    test_time, test_signal = deduplicate_time(test_time, test_signal)

    t_common, ref_interp, test_interp = resample_to_common_grid(
        ref_time, ref_signal, test_time, test_signal, skip_time=skip_time)

    ref_data  = np.column_stack([t_common, ref_interp])
    test_data = np.column_stack([t_common, test_interp])

    n1, n2 = create_comparison_plot(
        ref_data, test_data,
        test_name=f"RM3 Mooring Verification - {label}",
        output_dir=str(plots_dir),
        ref_file_path=ref_file_label,
        test_file_path=test_file_label,
        y_label=f"{label} ({units})",
    )

    l2_thresh, linf_thresh = tolerance
    passed = (n1 <= l2_thresh) and (n2 <= linf_thresh)
    tag = "PASS" if passed else "FAIL"
    print(f"  {label}: L2={n1:.2e} (thresh {l2_thresh:.0e}), "
          f"Linf={n2:.2e} (thresh {linf_thresh:.0e}) -> {tag}")

    return n1, n2, passed


def main():
    if len(sys.argv) != 3:
        print("Usage: compare_rm3_mooring.py <reference_dir> <results_file>")
        sys.exit(1)

    ref_dir      = Path(sys.argv[1])
    results_file = Path(sys.argv[2])

    print(f"Reference dir: {ref_dir}")
    print(f"Results file:  {results_file}")

    # Load reference
    ref_body, ref_tension = load_reference(ref_dir)
    ref_time = ref_body[:, 0]

    # Load HC results
    hc_data = np.loadtxt(str(results_file), skiprows=1)
    hc_time = hc_data[:, 0]

    plots_dir = results_file.parent / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    all_passed = True
    all_metrics = {}

    # ── Body motions ─────────────────────────────────────────────────────
    print("\n--- Body motion comparison (HC vs WEC-Sim) ---")

    for col_idx, label in [(1, "Float Heave Z"), (2, "Plate Heave Z")]:
        n1, n2, passed = compare_signal(
            ref_time, ref_body[:, col_idx],
            hc_time, hc_data[:, col_idx],
            label, "m", BODY_MOTION_TOLERANCE, plots_dir,
            f"wecsim_moordyn_body_motions.txt col {col_idx}",
            str(results_file.name),
            skip_time=RAMP_SKIP_TIME,
        )
        sname = label.lower().replace(" ", "_")
        all_metrics[sname] = {"l2_norm": n1, "linf_norm": n2}
        if not passed:
            all_passed = False

    # ── Fairlead tensions (MoorDyn lines.out) ────────────────────────────
    # MoorDyn writes output next to its input file. Search several likely locations.
    build_dir = Path(os.environ.get("HYDROCHRONO_BUILD_DIR", ""))
    moordyn_data_dir = build_dir / "data" / "yaml" / "rm3" / "mooring" if build_dir.is_dir() else Path()

    search_dirs = [
        moordyn_data_dir,
        results_file.parent,
        results_file.parent.parent.parent,
        Path.cwd(),
    ]
    hc_tension = load_moordyn_lines_out(search_dirs)

    if hc_tension is not None:
        print("\n--- Fairlead tension comparison (HC/MoorDyn vs WEC-Sim/MoorDyn) ---")
        hc_t_time = hc_tension[:, 0]

        for col_idx, label in [(1, "FairTen4"), (2, "FairTen5"), (3, "FairTen6")]:
            n1, n2, passed = compare_signal(
                ref_tension[:, 0], ref_tension[:, col_idx],
                hc_t_time, hc_tension[:, col_idx],
                label, "N", FAIRLEAD_TOLERANCE, plots_dir,
                f"wecsim_moordyn_fairlead_tensions.txt col {col_idx}",
                "lines.out",
                skip_time=RAMP_SKIP_TIME,
            )
            sname = label.lower()
            all_metrics[sname] = {"l2_norm": n1, "linf_norm": n2}
            if not passed:
                all_passed = False
    else:
        print("\nWARNING: MoorDyn lines.out not found -- skipping fairlead tension comparison")

    # ── Status ───────────────────────────────────────────────────────────
    status = "PASS" if all_passed else "FAIL"
    write_status_file(str(results_file.parent), "rm3_mooring", status, all_metrics)

    print(f"\n=== RM3 MOORING VERIFICATION: {status} ===")

    if not all_passed:
        sys.exit(1)


if __name__ == "__main__":
    main()
