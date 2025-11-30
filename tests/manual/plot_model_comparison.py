#!/usr/bin/env python3
"""
Model-Specific State-Space Comparison Figure

Creates a 3x3 grid comparison for a specific model (sphere or oswec):
  Row 1 (1x3): Kernel fit - RIRF vs State-Space fit for key DOF
  Row 2 (1x3): Time series - Convolution vs State-Space
  Row 3 (3x1 each): Performance scaling subplots

Usage:
    python plot_model_comparison.py sphere [--save]
    python plot_model_comparison.py oswec [--save]
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10


def load_model_data(model, output_dir):
    """Load all data files for a specific model."""
    data = {}
    
    # Kernel fit
    kernel_file = os.path.join(output_dir, f'{model}_kernel_fit.csv')
    if os.path.exists(kernel_file):
        data['kernel'] = pd.read_csv(kernel_file)
    
    # Time series
    conv_file = os.path.join(output_dir, f'{model}_decay_convolution.csv')
    ss_file = os.path.join(output_dir, f'{model}_decay_state_space.csv')
    if os.path.exists(conv_file) and os.path.exists(ss_file):
        data['conv'] = pd.read_csv(conv_file)
        data['ss'] = pd.read_csv(ss_file)
    
    # Performance scaling
    perf_file = os.path.join(output_dir, f'{model}_performance.csv')
    if os.path.exists(perf_file):
        data['performance'] = pd.read_csv(perf_file)
    
    return data


def plot_kernel_fit(ax, data, model):
    """Row 1: Kernel fit spanning all columns."""
    if 'kernel' not in data:
        ax.text(0.5, 0.5, f'No kernel data for {model}\n(Run {model}_decay_ss_comparison)', 
                ha='center', va='center', transform=ax.transAxes, fontsize=11)
        ax.set_title(f'{model.upper()}: Kernel Fit (Key DOF)')
        return
    
    df = data['kernel']
    ax.plot(df['t'], df['K_actual'], 'b-', label='RIRF Kernel', linewidth=2.5, alpha=0.9)
    ax.plot(df['t'], df['K_fit'], 'r--', label='State-Space Fit', linewidth=2.5, alpha=0.9)
    
    dof_name = 'Heave' if model == 'sphere' else 'Pitch'
    ax.set_xlabel('Time [s]', fontsize=11)
    ax.set_ylabel('K(t)', fontsize=11)
    ax.set_title(f'{model.upper()}: {dof_name} Kernel Fit (RIRF vs State-Space)', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, df['t'].max())
    
    # Calculate and show R²
    ss_res = np.sum((df['K_actual'] - df['K_fit'])**2)
    ss_tot = np.sum((df['K_actual'] - df['K_actual'].mean())**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    ax.text(0.98, 0.95, f'R² = {r2:.4f}', transform=ax.transAxes, ha='right', va='top',
            fontsize=11, bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))


def plot_time_series(ax, data, model):
    """Row 2: Time series comparison spanning all columns."""
    if 'conv' not in data or 'ss' not in data:
        ax.text(0.5, 0.5, f'No time series data for {model}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=11)
        return
    
    conv = data['conv']
    ss = data['ss']
    
    if model == 'sphere':
        y_conv = conv['heave_z']
        y_ss = ss['heave_z']
        ylabel = 'Heave [m]'
        title = 'Sphere Heave Decay'
        unit = 'm'
        scale = 1
    else:  # oswec
        y_conv = np.degrees(conv['flap_pitch_rad'])
        y_ss = np.degrees(ss['flap_pitch_rad'])
        ylabel = 'Pitch [deg]'
        title = 'OSWEC Flap Pitch Decay'
        unit = '°'
        scale = 1
    
    ax.plot(conv['time'], y_conv, 'b-', label='Convolution', linewidth=2, alpha=0.85)
    ax.plot(ss['time'], y_ss, 'r--', label='State-Space', linewidth=2, alpha=0.85)
    
    sim_time = conv['time'].max()
    ax.set_xlabel('Time [s]', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(f'{title} ({sim_time:.0f}s): Convolution vs State-Space', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    
    # Calculate and show max difference
    n = min(len(y_conv), len(y_ss))
    max_diff = np.max(np.abs(y_conv.values[:n] - y_ss.values[:n]))
    
    if model == 'sphere':
        diff_str = f'Max diff: {max_diff*100:.1f} cm'
    else:
        diff_str = f'Max diff: {max_diff:.2f}°'
    
    ax.text(0.02, 0.02, diff_str, transform=ax.transAxes, fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))


def plot_performance_total(ax, data, model):
    """Row 3, Col 1: Total computation time vs steps."""
    if 'performance' not in data:
        ax.text(0.5, 0.5, 'No performance data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Total Computation Time')
        return
    
    df = data['performance']
    steps = df['sim_steps'].values
    conv_time = df['conv_time_us'].values / 1000  # Convert to ms
    ss_time = df['ss_time_us'].values / 1000
    
    ax.loglog(steps, conv_time, 'b-o', label='Convolution', linewidth=2, markersize=7)
    ax.loglog(steps, ss_time, 'r-s', label='State-Space', linewidth=2, markersize=7)
    ax.set_xlabel('Simulation Steps', fontsize=10)
    ax.set_ylabel('Total Time [ms]', fontsize=10)
    ax.set_title('Total Computation Time', fontsize=11, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, alpha=0.3, which='both')


def plot_performance_per_step(ax, data, model):
    """Row 3, Col 2: Per-step cost showing O(N) vs O(1)."""
    if 'performance' not in data:
        ax.text(0.5, 0.5, 'No performance data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Per-Step Cost')
        return
    
    df = data['performance']
    steps = df['sim_steps'].values
    conv_per_step = df['conv_per_step_us'].values
    ss_per_step = df['ss_per_step_us'].values
    
    ax.semilogx(steps, conv_per_step, 'b-o', label='Convolution O(N)', linewidth=2, markersize=7)
    ax.semilogx(steps, ss_per_step, 'r-s', label='State-Space O(1)', linewidth=2, markersize=7)
    ax.set_xlabel('Simulation Steps', fontsize=10)
    ax.set_ylabel('Time per Step [μs]', fontsize=10)
    ax.set_title('Per-Step Cost', fontsize=11, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, alpha=0.3)
    
    # Annotate the scaling
    if len(steps) > 1:
        ratio = conv_per_step[-1] / conv_per_step[0]
        ax.annotate(f'Conv: {ratio:.0f}× slower\nat long sims', 
                   xy=(steps[-1], conv_per_step[-1]), 
                   xytext=(steps[len(steps)//2], conv_per_step[-1]*0.7),
                   fontsize=8, color='blue', ha='center',
                   arrowprops=dict(arrowstyle='->', color='blue', lw=1),
                   bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))


def plot_performance_speedup(ax, data, model):
    """Row 3, Col 3: Speedup factor bar chart."""
    if 'performance' not in data:
        ax.text(0.5, 0.5, 'No performance data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Speedup Factor')
        return
    
    df = data['performance']
    steps = df['sim_steps'].values
    speedup = df['speedup'].values
    
    colors = plt.cm.Greens(np.linspace(0.4, 0.8, len(steps)))
    bars = ax.bar(range(len(steps)), speedup, color=colors, edgecolor='black', linewidth=1.2)
    
    ax.set_xticks(range(len(steps)))
    ax.set_xticklabels([str(int(s)) for s in steps], fontsize=8, rotation=45, ha='right')
    ax.set_xlabel('Simulation Steps', fontsize=10)
    ax.set_ylabel('Speedup (Conv / SS)', fontsize=10)
    ax.set_title('State-Space Speedup', fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for bar, s in zip(bars, speedup):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{s:.0f}×', ha='center', va='bottom', fontsize=9, fontweight='bold')


def main():
    parser = argparse.ArgumentParser(description='Model-specific SS comparison figure')
    parser.add_argument('model', choices=['sphere', 'oswec'], help='Model to plot')
    parser.add_argument('--save', action='store_true', help='Save plot instead of displaying')
    parser.add_argument('--dir', default='output', help='Data directory')
    args = parser.parse_args()
    
    model = args.model
    output_dir = args.dir
    
    print(f"Creating comparison figure for: {model.upper()}")
    print(f"Data directory: {output_dir}")
    
    # Load data
    data = load_model_data(model, output_dir)
    print(f"  Kernel data: {'kernel' in data}")
    print(f"  Time series: {'conv' in data and 'ss' in data}")
    print(f"  Performance: {'performance' in data}")
    
    # Create figure with 3x3 grid
    fig = plt.figure(figsize=(14, 12))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 1], hspace=0.35, wspace=0.25)
    
    # Row 1: Kernel fit (spans all 3 columns)
    ax_kernel = fig.add_subplot(gs[0, :])
    plot_kernel_fit(ax_kernel, data, model)
    
    # Row 2: Time series (spans all 3 columns)
    ax_ts = fig.add_subplot(gs[1, :])
    plot_time_series(ax_ts, data, model)
    
    # Row 3: Performance subplots
    ax_perf1 = fig.add_subplot(gs[2, 0])
    ax_perf2 = fig.add_subplot(gs[2, 1])
    ax_perf3 = fig.add_subplot(gs[2, 2])
    
    plot_performance_total(ax_perf1, data, model)
    plot_performance_per_step(ax_perf2, data, model)
    plot_performance_speedup(ax_perf3, data, model)
    
    # Main title
    model_name = 'IEA Sphere' if model == 'sphere' else 'OSWEC'
    fig.suptitle(f'{model_name}: State-Space Radiation Damping Verification', 
                 fontsize=15, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if args.save:
        out_path = os.path.join(output_dir, f'{model}_ss_verification.png')
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"\nSaved: {out_path}")
    else:
        plt.show()
    
    plt.close()
    print("Done!")


if __name__ == '__main__':
    main()

