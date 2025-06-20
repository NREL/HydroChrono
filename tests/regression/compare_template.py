#!/usr/bin/env python3
"""
HydroChrono Regression Test Comparison Template

This template provides a standardized way to compare reference and test data
across different regression test cases. It generates professional comparison
plots with consistent formatting and comprehensive information panels.

Usage:
    python compare_template.py <reference_file> <test_file> [test_name] [y_label]

Example:
    python compare_template.py ref_data.txt test_data.txt "Sphere Decay Test" "Heave (m)"
"""

import sys
import os
from pathlib import Path
from datetime import datetime
import platform

# Check for required packages
try:
    import numpy as np
except ImportError:
    print("Error: numpy is required but not installed. Please install it with: pip install numpy")
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: matplotlib is required but not installed. Please install it with: pip install matplotlib")
    sys.exit(1)

# Layout configuration for text boxes and panels - using absolute positioning
LAYOUT = {
    'figure': {
        'figsize': (12, 9),
        'facecolor': 'white'
    },
    'fonts': {
        'title': 13,
        'heading': 11,
        'body': 11,
        'small': 8
    },
    'panels': {
        'test_info': {
            'pos': (0.02, 0.82, 0.8, 0.12),
            'font_size': 'body',
            'style': {
                'facecolor': '#f8f9fa',
                'edgecolor': '#e9ecef',
                'text_color': '#212529'
            }
        },
        'system_info': {
            'pos': (0.85, 0.82, 0.22, 0.12),
            'font_size': 'body',
            'style': {
                'facecolor': '#f8f9fa',
                'edgecolor': '#e9ecef',
                'text_color': '#212529'
            }
        },
        'simulation_summary': {
            'pos': (0.85, 0.63, 0.22, 0.15),
            'font_size': 'body',
            'style': {
                'facecolor': '#f8f9fa',
                'edgecolor': '#e9ecef',
                'text_color': '#212529'
            }
        },
        'data_stats': {
            'pos': (0.85, 0.18, 0.22, 0.35),
            'font_size': 'body',
            'style': {
                'facecolor': '#f8f9fa',
                'edgecolor': '#e9ecef',
                'text_color': '#212529'
            }
        },
        'error_metrics': {
            'pos': (0.85, 0.07, 0.22, 0.15),
            'font_size': 'body',
            'style': {
                'facecolor': '#f8f9fa',
                'edgecolor': '#e9ecef',
                'text_color': '#dc3545'
            }
        }
    },
    'plots': {
        'main_plot': {
            'pos': (0.05, 0.45, 0.78, 0.28)
        },
        'error_plot': {
            'pos': (0.05, 0.05, 0.78, 0.28)
        }
    }
}

def format_path(path):
    """Format file paths for display by making them relative to current directory"""
    if path is None:
        return "Not specified"
    try:
        rel_path = os.path.relpath(path, os.getcwd())
        return rel_path if len(rel_path) < len(path) else path
    except (ValueError, OSError):
        return str(path)

def get_cmake_cache_path():
    """Return the canonical path to CMakeCache.txt in the build directory."""
    # Try from current working directory first (for CTest running from results dir)
    cwd_path = os.path.join(os.getcwd(), '../../../../../CMakeCache.txt')
    if os.path.exists(cwd_path):
        return cwd_path
    
    # Fallback to script location relative to project root
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    cmake_cache_path = os.path.join(project_root, 'build', 'CMakeCache.txt')
    return cmake_cache_path if os.path.exists(cmake_cache_path) else None

def get_hydrochrono_version():
    """Get HydroChrono version from CMakeCache.txt"""
    try:
        cmake_cache_path = get_cmake_cache_path()
        if cmake_cache_path:
            with open(cmake_cache_path, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('CMAKE_PROJECT_VERSION:STATIC='):
                        return line.split('=')[1].strip()
        return os.environ.get('HYDROCHRONO_VERSION', 'Unknown')
    except (OSError, IOError, UnicodeDecodeError):
        return os.environ.get('HYDROCHRONO_VERSION', 'Unknown')

def get_chrono_version():
    """Get Chrono version from Chrono CMakeLists.txt"""
    try:
        cmake_cache_path = get_cmake_cache_path()
        if not cmake_cache_path:
            return os.environ.get('CHRONO_VERSION', 'Unknown')
        
        # Get Chrono_DIR from CMakeCache.txt
        chrono_dir = None
        with open(cmake_cache_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('Chrono_DIR:UNINITIALIZED='):
                    chrono_dir = line.split('=')[1].strip()
                    break
        
        if chrono_dir:
            # Navigate to Chrono root directory
            chrono_root = os.path.dirname(os.path.dirname(chrono_dir))
            chrono_cmakelists = os.path.join(chrono_root, 'CMakeLists.txt')
            
            if os.path.exists(chrono_cmakelists):
                with open(chrono_cmakelists, 'r', encoding='utf-8') as f:
                    major = minor = patch = "0"
                    for line in f:
                        line = line.strip()
                        if line.startswith('set(CHRONO_VERSION_MAJOR'):
                            # Extract number from: set(CHRONO_VERSION_MAJOR 9)
                            parts = line.split()
                            if len(parts) >= 2:
                                major = parts[1].rstrip(')')
                        elif line.startswith('set(CHRONO_VERSION_MINOR'):
                            parts = line.split()
                            if len(parts) >= 2:
                                minor = parts[1].rstrip(')')
                        elif line.startswith('set(CHRONO_VERSION_PATCH'):
                            parts = line.split()
                            if len(parts) >= 2:
                                patch = parts[1].rstrip(')')
                    
                    if major != "0" or minor != "0" or patch != "0":
                        return f"{major}.{minor}.{patch}"
        
        return os.environ.get('CHRONO_VERSION', 'Unknown')
    except (OSError, IOError, UnicodeDecodeError):
        return os.environ.get('CHRONO_VERSION', 'Unknown')

def create_text_panel(fig, panel_config, content):
    """Create a text panel with given configuration using absolute positioning"""
    left, bottom, width, height = panel_config['pos']
    ax = fig.add_axes([left, bottom, width, height])
    ax.axis('off')
    style = panel_config['style']
    font_size = LAYOUT['fonts'][panel_config['font_size']]
    
    return ax.text(0.05, 0.95, content,
                  transform=ax.transAxes,
                  verticalalignment='top',
                  bbox=dict(boxstyle='round,pad=0.6',
                           facecolor=style['facecolor'],
                           edgecolor=style['edgecolor'],
                           linewidth=1.5,
                           alpha=0.95),
                  fontsize=font_size,
                  family='monospace',
                  fontweight='normal',
                  color=style['text_color'])

def apply_modern_style(ax):
    """Apply consistent modern styling to plot axes"""
    ax.grid(True, alpha=0.2, color='#6c757d', linewidth=0.5)
    ax.tick_params(labelsize=9, colors='#495057')
    for spine in ['top', 'right', 'left', 'bottom']:
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_color('#dee2e6')
        ax.spines[spine].set_linewidth(1.0)
    ax.set_facecolor('#ffffff')

def find_executable(test_dir, executable_patterns):
    """
    Find executable in test directory or its parent
    
    Args:
        test_dir: Directory to search in
        executable_patterns: List of patterns to search for (e.g., ["sphere_decay_test", "rm3_test"])
    
    Returns:
        Path to executable if found, None otherwise
    """
    search_dirs = [test_dir, test_dir.parent]

    try:
        for s_dir in search_dirs:
            for pattern in executable_patterns:
                # Look for common executable names
                possible_names = [pattern, f"{pattern}.exe", f"{pattern}.out"]
                
                for name in possible_names:
                    exe_file = s_dir / name
                    if exe_file.exists():
                        return exe_file
                
                # If not found, look for any executable with the pattern in the name
                for exe_file in s_dir.glob("*"):
                    if exe_file.is_file() and pattern in exe_file.name:
                        # Check if it's executable (Unix) or has executable extension (Windows)
                        if (os.access(exe_file, os.X_OK) or 
                            exe_file.suffix in ['.exe', '.out', '.app']):
                            return exe_file
    except Exception as e:
        print(f"Warning: Could not find executable: {e}")
    
    return None

def create_comparison_plot(ref_data, test_data, test_name, output_dir, 
                          ref_file_path=None, test_file_path=None, executable_path=None,
                          y_label="Value", executable_patterns=None):
    """
    Create comparison plot between reference and test data
    
    Args:
        ref_data: Reference data array
        test_data: Test data array  
        test_name: Name of the test case
        output_dir: Directory to save the plot
        ref_file_path: Path to reference file (for transparency)
        test_file_path: Path to test file (for transparency)
        executable_path: Path to the executable that generated the test data
        y_label: Label for the y-axis (e.g., "Heave (m)", "Surge (m)")
        executable_patterns: List of patterns to search for executable if not provided
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Calculate error metrics
    nval = test_data.shape[0]
    x = np.linspace(test_data[0, 0], test_data[nval-1, 0], nval)
    y1 = np.interp(x, ref_data[:,0], ref_data[:,1])
    y2 = np.interp(x, test_data[:,0], test_data[:,1])
    yd = y1 - y2
    n1 = np.linalg.norm(yd)/nval  # L2 norm
    n2 = np.linalg.norm(yd, np.inf)  # L-infinity norm
    
    # Set up the figure with configuration
    fig_cfg = LAYOUT['figure']
    fig = plt.figure(figsize=fig_cfg['figsize'], facecolor=fig_cfg['facecolor'])
    
    # Extract model name from executable path, falling back to test file stem
    model_name = "Unknown Model"
    if executable_path:
        exe_name = os.path.basename(executable_path)
        # Remove common executable extensions
        for ext in ['.exe', '.out', '.app']:
            if exe_name.endswith(ext):
                model_name = exe_name[:-len(ext)]
                break
        else:
            model_name = exe_name
    elif test_file_path:
        model_name = Path(test_file_path).stem
    
    # Create Test Information panel
    info_content = (
        f"Test Information\n\n"
        f"Model/Executable: {model_name}\n"
        f"Reference File: {format_path(ref_file_path)}\n"
        f"Simulation File: {format_path(test_file_path)}\n"
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    )
    create_text_panel(fig, LAYOUT['panels']['test_info'], info_content)
    
    # Create System Information panel
    platform_name = {"Windows": "Windows", "Darwin": "macOS"}.get(platform.system(), "Linux")
    python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    hydrochrono_version = get_hydrochrono_version()
    chrono_version = get_chrono_version()
    
    sysinfo_content = (
        f"System Information\n\n"
        f"Platform: {platform_name}\n"
        f"Python: {python_version}\n"
        f"HydroChrono: {hydrochrono_version}\n"
        f"Chrono: {chrono_version}"
    )
    create_text_panel(fig, LAYOUT['panels']['system_info'], sysinfo_content)
    
    # Create main comparison plot
    plot_cfg = LAYOUT['plots']['main_plot']
    ax1 = fig.add_axes(plot_cfg['pos'])
    ax1.plot(ref_data[:,0], ref_data[:,1], color='#007bff', linewidth=2.5, label='Reference', alpha=0.9)
    ax1.plot(test_data[:,0], test_data[:,1], color='#dc3545', linewidth=2.5, label='Simulation', alpha=0.9, linestyle='--')
    ax1.set_xlabel('Time (s)', fontsize=11, color='#495057', fontweight='500')
    ax1.set_ylabel(y_label, fontsize=11, color='#495057', fontweight='500')
    ax1.set_title(f'{test_name} - Reference vs Simulation Comparison', fontsize=13, fontweight='bold', color='#212529', pad=25)
    ax1.legend(fontsize=10, framealpha=0.9)
    apply_modern_style(ax1)
    
    # Create error plot
    plot_cfg = LAYOUT['plots']['error_plot']
    ax2 = fig.add_axes(plot_cfg['pos'])
    ax2.plot(x, yd, color='#dc3545', linewidth=2, label='Error (Ref - Sim)', alpha=0.8)
    ax2.axhline(y=0, color='#6c757d', linestyle='-', alpha=0.4, linewidth=1)
    ax2.set_xlabel('Time (s)', fontsize=11, color='#495057', fontweight='500')
    ax2.set_ylabel('Error (m)', fontsize=11, color='#495057', fontweight='500')
    ax2.set_title('Error Analysis', fontsize=13, fontweight='bold', color='#212529', pad=15)
    ax2.legend(fontsize=10, framealpha=0.9)
    apply_modern_style(ax2)
    
    # Create Error Metrics panel
    content_text = (
        f"Error Metrics\n\n"
        f"L₂ Norm: {n1:.2e}\n"
        f"L∞ Norm: {n2:.2e}\n"
        f"Max Error: {np.max(np.abs(yd)):.2e}\n"
        f"Mean Error: {np.mean(yd):.2e}"
    )
    create_text_panel(fig, LAYOUT['panels']['error_metrics'], content_text)
    
    # Create Data Statistics panel
    stats_content = (
        f"Data Statistics\n\n"
        f"Reference:\n"
        f"  Mean: {np.mean(ref_data[:,1]):.4f}m\n"
        f"  Std: {np.std(ref_data[:,1]):.4f}m\n"
        f"  Range: [{np.min(ref_data[:,1]):.3f}, {np.max(ref_data[:,1]):.3f}]m\n\n"
        f"Simulation:\n"
        f"  Mean: {np.mean(test_data[:,1]):.4f}m\n"
        f"  Std: {np.std(test_data[:,1]):.4f}m\n"
        f"  Range: [{np.min(test_data[:,1]):.3f}, {np.max(test_data[:,1]):.3f}]m\n\n"
        f"Correlation: {np.corrcoef(ref_data[:,1], test_data[:,1])[0,1]:.6f}"
    )
    create_text_panel(fig, LAYOUT['panels']['data_stats'], stats_content)
    
    # Create Simulation Summary panel
    ref_duration = ref_data[-1, 0] - ref_data[0, 0]
    test_duration = test_data[-1, 0] - test_data[0, 0]
    ref_dt = ref_data[1, 0] - ref_data[0, 0] if len(ref_data) > 1 else 0
    test_dt = test_data[1, 0] - test_data[0, 0] if len(test_data) > 1 else 0
    
    summary_content = (
        f"Simulation Summary\n\n"
        f"Reference:\n"
        f"  Duration: {ref_duration:.1f}s\n"
        f"  Timestep: {ref_dt:.3f}s\n"
        f"  Points: {len(ref_data)}\n\n"
        f"Simulation:\n"
        f"  Duration: {test_duration:.1f}s\n"
        f"  Timestep: {test_dt:.3f}s\n"
        f"  Points: {len(test_data)}"
    )
    create_text_panel(fig, LAYOUT['panels']['simulation_summary'], summary_content)
    
    # Save plot
    plot_filename = os.path.join(output_dir, f'{test_name}_comparison.png')
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    print(f"Plot saved: {plot_filename}")
    
    return n1, n2

def run_comparison(ref_file, test_file, test_name=None, y_label="Value", 
                  executable_patterns=None, pass_criteria=None):
    """
    Run a complete comparison between reference and test data
    
    Args:
        ref_file: Path to reference data file
        test_file: Path to test data file
        test_name: Name of the test (defaults to test file stem)
        y_label: Label for y-axis
        executable_patterns: List of patterns to search for executable
        pass_criteria: Tuple of (l2_threshold, linf_threshold) for pass/fail
    
    Returns:
        Tuple of (l2_norm, linf_norm)
    """
    print(f"Comparing: {ref_file} vs {test_file}")

    # Load data with error handling
    try:
        refData = np.loadtxt(ref_file, skiprows=1)
        testData = np.loadtxt(test_file, skiprows=1)
    except (OSError, IOError, ValueError) as e:
        print(f"Error loading data files: {e}")
        sys.exit(1)
    
    # Validate data
    if refData.size == 0 or testData.size == 0:
        print("Error: One or both data files are empty")
        sys.exit(1)
    
    if refData.shape[1] < 2 or testData.shape[1] < 2:
        print("Error: Data files must have at least 2 columns (time and value)")
        sys.exit(1)

    print(f"Reference data shape: {refData.shape}")
    print(f"Test data shape: {testData.shape}")

    # Determine test name from file path if not provided
    if test_name is None:
        test_name = Path(test_file).stem
    
    # Create plots directory in the same location as the test file
    test_file_path = Path(test_file)
    plots_dir = test_file_path.parent / "plots"
    
    # Find the executable path
    executable_path = None
    if executable_patterns:
        executable_path = find_executable(test_file_path.parent, executable_patterns)
    
    # Generate comparison plot
    def rel_to_root(path):
        try:
            project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
            return os.path.relpath(path, project_root)
        except Exception:
            return str(path)
    
    try:
        n1, n2 = create_comparison_plot(
            refData, testData, test_name, plots_dir, 
            ref_file_path=rel_to_root(ref_file), 
            test_file_path=rel_to_root(test_file),
            executable_path=rel_to_root(str(executable_path)) if executable_path else None,
            y_label=y_label,
            executable_patterns=executable_patterns
        )
    except Exception as e:
        print(f"Error creating comparison plot: {e}")
        sys.exit(1)
    
    # Check pass/fail criteria if provided
    if pass_criteria:
        l2_threshold, linf_threshold = pass_criteria
        if (n1 > l2_threshold or n2 > linf_threshold):
            print(f"TEST FAILED - L2 Norm: {n1:.2e}, L-infinity Norm: {n2:.2e}")
            return n1, n2, False
        else:
            print(f"TEST PASSED - L2 Norm: {n1:.2e}, L-infinity Norm: {n2:.2e}")
            return n1, n2, True
    
    return n1, n2

if __name__ == '__main__':
    """
    Template usage example - customize this section for your specific test case
    """
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    ref_file = sys.argv[1]
    test_file = sys.argv[2]
    test_name = sys.argv[3] if len(sys.argv) > 3 else None
    y_label = sys.argv[4] if len(sys.argv) > 4 else "Value"
    
    # Customize these for your specific test case
    executable_patterns = ["test_executable"]  # Add patterns to search for
    pass_criteria = (1e-4, 0.02)  # (L2 threshold, L-infinity threshold)
    
    run_comparison(ref_file, test_file, test_name, y_label, 
                  executable_patterns, pass_criteria) 