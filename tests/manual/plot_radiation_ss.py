#!/usr/bin/env python3
"""
Plot results from radiation_ss_visual_check.

Usage:
    python plot_radiation_ss.py [data_dir]

    If data_dir is not specified, looks for CSV files in ./output/ subdirectory.

Requires: matplotlib, pandas (numpy comes with pandas)
    pip install matplotlib pandas
"""

import sys
from pathlib import Path

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Install with: pip install matplotlib pandas")
    sys.exit(1)


def plot_constant_velocity(data_dir: Path, save: bool = False):
    """Plot constant velocity response with model vs exact comparison."""
    csv_file = data_dir / "radiation_ss_constant_v.csv"
    if not csv_file.exists():
        print(f"  Skipping: {csv_file} not found")
        return

    df = pd.read_csv(csv_file)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Constant Velocity Response (α=0.5, H=2.0, v₀=1.0)", fontsize=14, fontweight='bold')

    # Force comparison
    ax = axes[0, 0]
    ax.plot(df["t"], df["f_model"], "b-", linewidth=2, label="Model")
    ax.plot(df["t"], df["f_exact"], "r--", linewidth=2, label="Exact")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radiation Force")
    ax.set_title("Force: Model vs Exact")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # State comparison
    ax = axes[0, 1]
    ax.plot(df["t"], df["z_model"], "b-", linewidth=2, label="Model")
    ax.plot(df["t"], df["z_exact"], "r--", linewidth=2, label="Exact")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Internal State z")
    ax.set_title("State: Model vs Exact")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Error plot
    ax = axes[1, 0]
    ax.semilogy(df["t"][1:], df["error_pct"][1:], "g-", linewidth=1.5)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Relative Error [%]")
    ax.set_title("Relative Force Error")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=1e-6)

    # Phase portrait (z vs f)
    ax = axes[1, 1]
    ax.plot(df["z_model"], df["f_model"], "b-", linewidth=2, label="Model trajectory")
    ax.plot(df["z_exact"], df["f_exact"], "r--", linewidth=2, label="Exact trajectory")
    ax.scatter([df["z_model"].iloc[-1]], [df["f_model"].iloc[-1]], 
               c='b', s=100, marker='o', zorder=5, label="Final state (model)")
    ax.set_xlabel("Internal State z")
    ax.set_ylabel("Radiation Force")
    ax.set_title("Phase Portrait")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    
    if save:
        out_file = data_dir / "plot_constant_velocity.png"
        plt.savefig(out_file, dpi=150)
        print(f"  Saved: {out_file}")
    
    plt.show()


def plot_step_velocity(data_dir: Path, save: bool = False):
    """Plot step velocity response."""
    csv_file = data_dir / "radiation_ss_step_v.csv"
    if not csv_file.exists():
        print(f"  Skipping: {csv_file} not found")
        return

    df = pd.read_csv(csv_file)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle("Step Velocity Response (v: 1→0 at t=5s, α=1.0, H=5.0)", 
                 fontsize=14, fontweight='bold')

    # Velocity input
    ax = axes[0]
    ax.plot(df["t"], df["v"], "k-", linewidth=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity v")
    ax.set_title("Input Velocity")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 1.2)

    # Internal state
    ax = axes[1]
    ax.plot(df["t"], df["z_model"], "b-", linewidth=2)
    ax.axvline(x=5.0, color='gray', linestyle='--', alpha=0.7, label='Step time')
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Internal State z")
    ax.set_title("Internal State Response")
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Force output
    ax = axes[2]
    ax.plot(df["t"], df["f_model"], "r-", linewidth=2)
    ax.axvline(x=5.0, color='gray', linestyle='--', alpha=0.7, label='Step time')
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radiation Force")
    ax.set_title("Force Response")
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    
    if save:
        out_file = data_dir / "plot_step_velocity.png"
        plt.savefig(out_file, dpi=150)
        print(f"  Saved: {out_file}")
    
    plt.show()


def plot_multi_mode(data_dir: Path, save: bool = False):
    """Plot multiple mode superposition."""
    csv_file = data_dir / "radiation_ss_multi_mode.csv"
    if not csv_file.exists():
        print(f"  Skipping: {csv_file} not found")
        return

    df = pd.read_csv(csv_file)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle("Multiple Modes Superposition", fontsize=14, fontweight='bold')

    # Individual mode states
    ax = axes[0]
    ax.plot(df["t"], df["z0_model"], "b-", linewidth=2, label="Mode 0 (slow, τ=2s)")
    ax.plot(df["t"], df["z1_model"], "orange", linewidth=2, label="Mode 1 (fast, τ=0.5s)")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Internal State z")
    ax.set_title("Mode States")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Individual mode forces
    ax = axes[1]
    ax.plot(df["t"], df["f_mode0"], "b-", linewidth=2, label="Mode 0 (H=3)")
    ax.plot(df["t"], df["f_mode1"], "orange", linewidth=2, label="Mode 1 (H=1)")
    ax.plot(df["t"], df["f0_exact"], "b--", linewidth=1, alpha=0.7)
    ax.plot(df["t"], df["f1_exact"], "--", color='orange', linewidth=1, alpha=0.7)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Force Contribution")
    ax.set_title("Mode Force Contributions")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Total force
    ax = axes[2]
    ax.plot(df["t"], df["f_total"], "g-", linewidth=2, label="Model total")
    ax.plot(df["t"], df["f_total_exact"], "r--", linewidth=2, label="Exact total")
    ax.axhline(y=6.0+0.5, color='gray', linestyle=':', alpha=0.5, label='Steady state (6.5)')
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Total Radiation Force")
    ax.set_title("Total Force: Model vs Exact")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    
    if save:
        out_file = data_dir / "plot_multi_mode.png"
        plt.savefig(out_file, dpi=150)
        print(f"  Saved: {out_file}")
    
    plt.show()


def plot_dt_sensitivity(data_dir: Path, save: bool = False):
    """Plot time step sensitivity analysis."""
    csv_file = data_dir / "radiation_ss_dt_sensitivity.csv"
    if not csv_file.exists():
        print(f"  Skipping: {csv_file} not found")
        return

    df = pd.read_csv(csv_file)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle("Time Step Sensitivity (α=2.0, H=5.0)", fontsize=14, fontweight='bold')

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    labels = ['dt=0.01s', 'dt=0.1s', 'dt=0.5s', 'dt=1.0s']
    cols = ['f_dt001', 'f_dt01', 'f_dt05', 'f_dt1']

    # Force comparison
    ax = axes[0]
    for col, label, color in zip(cols, labels, colors):
        ax.plot(df["t"], df[col], color=color, linewidth=1.5, label=label)
    ax.plot(df["t"], df["f_exact"], "k--", linewidth=2, label="Exact", alpha=0.7)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radiation Force")
    ax.set_title("Force Response with Different Time Steps")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Error comparison
    ax = axes[1]
    for col, label, color in zip(cols, labels, colors):
        error = np.abs(df[col] - df["f_exact"]) / (np.abs(df["f_exact"]) + 1e-10) * 100
        ax.semilogy(df["t"][1:], error[1:], color=color, linewidth=1.5, label=label)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Relative Error [%]")
    ax.set_title("Relative Error vs Time Step")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(1e-8, 100)

    plt.tight_layout()
    
    if save:
        out_file = data_dir / "plot_dt_sensitivity.png"
        plt.savefig(out_file, dpi=150)
        print(f"  Saved: {out_file}")
    
    plt.show()


def main():
    # Parse command line arguments
    save = "--save" in sys.argv or "-s" in sys.argv
    
    # Filter out flags to find the data directory argument
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    data_dir = Path(args[0]) if args else Path("output")
    
    print("=" * 60)
    print("RadiationStateSpaceModel Visual Verification Plots")
    print("=" * 60)
    print(f"Data directory: {data_dir.absolute()}")
    print(f"Save plots: {save}")
    print()

    # Check if any CSV files exist
    csv_files = list(data_dir.glob("radiation_ss_*.csv"))
    if not csv_files:
        print(f"No CSV files found in {data_dir}")
        print("Run radiation_ss_visual_check first to generate data.")
        return 1

    print(f"Found {len(csv_files)} data file(s):")
    for f in csv_files:
        print(f"  - {f.name}")
    print()

    # Plot each scenario
    print("Plotting constant velocity response...")
    plot_constant_velocity(data_dir, save)
    
    print("Plotting step velocity response...")
    plot_step_velocity(data_dir, save)
    
    print("Plotting multiple mode response...")
    plot_multi_mode(data_dir, save)
    
    print("Plotting time step sensitivity...")
    plot_dt_sensitivity(data_dir, save)

    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(main())

