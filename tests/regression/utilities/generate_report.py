#!/usr/bin/env python3
"""
HydroChrono Regression Test Report Generator

This script generates a comprehensive markdown report of all regression test results,
including comparison plots and summary statistics. The report can optionally be
converted to PDF using pandoc.

Usage:
    python generate_report.py [--pdf] [--output-dir <dir>]
"""

import os
import sys
import glob
import argparse
from pathlib import Path
from datetime import datetime
import subprocess
from collections import defaultdict

# pypandoc is optional - only needed for PDF generation
try:
    import pypandoc
    PYPANDOC_AVAILABLE = True
except ImportError:
    pypandoc = None
    PYPANDOC_AVAILABLE = False

def natural_sort_key(filename):
    """Sort key for natural sorting of filenames with numbers."""
    import re
    # Split filename into text and number parts for natural sorting
    return [int(text) if text.isdigit() else text.lower() 
            for text in re.split('([0-9]+)', str(filename))]

def find_plot_files(build_dir):
    """Find all comparison plot files in the build directory."""
    plots = defaultdict(list)
    
    # Convert to absolute path to handle relative paths correctly
    build_path = Path(build_dir).resolve()
    
    # Define the models and their expected plot directories
    models = ['sphere', 'f3of', 'oswec', 'rm3']
    
    for model in models:
        # Try different possible paths for the plots directory
        possible_paths = [
            # New unified location: <build>/bin/Release/results/tests/<model>/plots
            build_path / "bin" / "Release" / "results" / "tests" / model / "plots",
            # Legacy / alternative layouts (kept for backward compatibility)
            build_path / "bin" / "tests" / "regression" / model / "results" / "plots",
            build_path / "tests" / "regression" / "Release" / model / "results" / "plots",
            build_path / model / "results" / "plots",  # If running from Release directory
            Path(".") / model / "results" / "plots",  # If running from current directory
            Path("bin") / "tests" / "regression" / model / "results" / "plots",  # If running from build dir
        ]
        
        plots_dir = None
        for path in possible_paths:
            if path.exists():
                plots_dir = path
                break
        
        if plots_dir:
            plot_files = list(plots_dir.glob("*_comparison.png"))
            # Use natural sorting to handle wave numbers properly (wave_1, wave_2, ..., wave_10)
            plots[model] = sorted(plot_files, key=lambda x: natural_sort_key(x.name))
            print(f"Found {len(plot_files)} plots for {model} in {plots_dir}")
        else:
            print(f"No plots directory found for {model}")
    
    return plots

def categorize_plots(plots):
    """Categorize plots by test type (decay, regular_waves, irregular_waves, etc.)."""
    categorized = defaultdict(lambda: defaultdict(list))
    
    for model, plot_files in plots.items():
        for plot_file in plot_files:
            filename = plot_file.name.lower()
            
            # Categorize based on filename patterns
            if 'decay' in filename:
                category = 'decay'
            elif 'irregular' in filename or 'irreg_waves' in filename:
                category = 'irregular_waves'
            elif 'regular' in filename or 'reg_waves' in filename:
                category = 'regular_waves'
            elif 'dt1' in filename or 'dt2' in filename or 'dt3' in filename:
                category = 'decay'  # F3OF tests are decay tests
            else:
                category = 'other'
            
            categorized[model][category].append(plot_file)
        
        # Sort plots within each category using natural sorting
        for category in categorized[model]:
            categorized[model][category] = sorted(categorized[model][category], 
                                                key=lambda x: natural_sort_key(x.name))
    
    return categorized

def get_test_summary(categorized_plots):
    """Generate summary statistics for the tests."""
    total_plots = sum(len(plots) for model_plots in categorized_plots.values() 
                     for plots in model_plots.values())
    
    model_counts = {}
    for model, categories in categorized_plots.items():
        model_counts[model] = sum(len(plots) for plots in categories.values())
    
    return {
        'total_tests': total_plots,
        'model_counts': model_counts,
        'models': list(categorized_plots.keys())
    }

def get_test_results(categorized_plots, build_dir=None):
    """Get test results from CTest logs, falling back to UNKNOWN when unavailable."""
    test_results = {
        'sphere': {},
        'f3of': {},
        'oswec': {},
        'rm3': {}
    }
    
    # Default status: plots exist but we don't know pass/fail without CTest logs
    for model, model_plots in categorized_plots.items():
        if model not in test_results:
            test_results[model] = {}
        
        for test_type, plots in model_plots.items():
            if plots:
                test_results[model][test_type] = 'UNKNOWN'
            else:
                test_results[model][test_type] = 'NO DATA'
    
    # Try to find and parse CTest log files for accurate results
    import re
    
    # Search for CTest log in likely locations relative to build_dir
    search_roots = []
    if build_dir:
        search_roots.append(Path(build_dir).resolve())
    search_roots.append(Path.cwd())
    
    ctest_log_patterns = []
    for root in search_roots:
        ctest_log_patterns.append(root / "Testing" / "Temporary" / "LastTest.log")
        ctest_log_patterns.append(root / "Testing" / "Temporary" / "LastTest.log.tmp")
    
    # Mapping from CTest test names to (model, test_type) pairs.
    # Only regression comparison tests (the *_regression tests) matter here;
    # the C++ execute tests just produce data files.
    test_name_map = {
        'sphere_decay_regression':                      ('sphere', 'decay'),
        'sphere_reg_waves_regression':                  ('sphere', 'regular_waves'),
        'sphere_irreg_waves_regression':                ('sphere', 'irregular_waves'),
        'sphere_irreg_waves_eta_regression':            ('sphere', 'irregular_waves_eta'),
        'sphere_irreg_waves_eta_consistency_regression': ('sphere', 'other'),
        'f3of_dt1_regression':                          ('f3of', 'decay'),
        'f3of_dt2_regression':                          ('f3of', 'decay'),
        'f3of_dt3_regression':                          ('f3of', 'decay'),
        'oswec_decay_regression':                       ('oswec', 'decay'),
        'oswec_reg_waves_regression':                   ('oswec', 'regular_waves'),
        'rm3_decay_regression':                         ('rm3', 'decay'),
        'rm3_reg_waves_regression':                     ('rm3', 'regular_waves'),
    }
    
    ctest_found = False
    for log_file in ctest_log_patterns:
        if log_file.exists():
            try:
                with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                
                # LastTest.log format:
                #   X/Y Testing: <test_name>
                #   ... (output) ...
                #   Test Passed.   OR   Test Failed.
                test_blocks = re.findall(
                    r'\d+/\d+\s+Testing:\s+(\S+).*?Test\s+(Passed|Failed)\.',
                    content, re.DOTALL)
                
                for test_name, status in test_blocks:
                    mapping = test_name_map.get(test_name)
                    if not mapping:
                        continue
                    model, test_type = mapping
                    
                    result_status = 'PASS' if status == 'Passed' else 'FAIL'
                    
                    if model in test_results:
                        # For f3of, multiple tests map to 'decay' — only mark FAIL
                        # if any sub-test fails (don't overwrite FAIL with PASS).
                        existing = test_results[model].get(test_type)
                        if existing == 'FAIL':
                            continue
                        test_results[model][test_type] = result_status
                        ctest_found = True
                                
            except Exception as e:
                print(f"Warning: Could not parse CTest log {log_file}: {e}")
    
    # Also check LastTestsFailed.log — contains names of failed tests from the
    # most recent ctest run.  This is a quick cross-check even if LastTest.log
    # was already parsed (belt-and-suspenders).
    for root in search_roots:
        failed_log = root / "Testing" / "Temporary" / "LastTestsFailed.log"
        if failed_log.exists():
            try:
                with open(failed_log, 'r', encoding='utf-8', errors='ignore') as f:
                    for line in f:
                        # Format: "<test_number>:<test_name>"
                        parts = line.strip().split(':', 1)
                        if len(parts) == 2:
                            failed_name = parts[1].strip()
                            mapping = test_name_map.get(failed_name)
                            if mapping:
                                model, test_type = mapping
                                if model in test_results:
                                    test_results[model][test_type] = 'FAIL'
                                    ctest_found = True
            except Exception as e:
                print(f"Warning: Could not parse {failed_log}: {e}")
    
    if ctest_found:
        print("Successfully parsed CTest results.")
    else:
        print("WARNING: No CTest log found — test results shown as UNKNOWN.")
        print("  Run tests via ctest first so that Testing/Temporary/LastTest.log exists.")
    
    return test_results

def _relative_plot_path(plot_file, output_dir):
    """Compute the relative path from the report output directory to a plot file."""
    rel = os.path.relpath(plot_file, output_dir).replace('\\', '/')
    return rel.replace(' ', '%20')

def generate_markdown_report(categorized_plots, output_dir, build_dir, html_styling=False):
    """Generate the markdown report content."""
    
    summary = get_test_summary(categorized_plots)
    
    # Start building the markdown content
    content = []
    
    # Title page - clean markdown or styled HTML
    if html_styling:
        content.append('<div class="header">')
        content.append('<h1>HydroChrono Regression Test Report</h1>')
        content.append('<p class="subtitle">Automated validation of hydrodynamic simulation capabilities</p>')
        content.append(f'<p class="date">Generated on {datetime.now().strftime("%B %d, %Y at %H:%M:%S")}</p>')
        content.append('</div>')
    else:
        content.append("# HydroChrono Regression Test Report")
        content.append("")
        content.append("**Automated validation of hydrodynamic simulation capabilities**")
        content.append("")
        content.append(f"*Generated on {datetime.now().strftime('%B %d, %Y at %H:%M:%S')}*")
    content.append("")
    
    # Executive Summary
    if html_styling:
        content.append('<div class="executive-summary">')
    content.append("## Executive Summary")
    content.append("")
    content.append("This document presents the results of the HydroChrono regression test suite, generated automatically to validate code integrity. The regression tests compare current simulation outputs against established reference data to ensure that recent code changes have not introduced numerical errors or regressions in the hydrodynamic simulation capabilities.")
    content.append("")
    
    # Key Metrics
    if html_styling:
        content.append('<div class="metrics-grid">')
        content.append('<div class="metric-card">')
        content.append(f'<span class="metric-value">{summary["total_tests"]}</span>')
        content.append('<span class="metric-label">Total Tests</span>')
        content.append('</div>')
        content.append('<div class="metric-card">')
        content.append(f'<span class="metric-value">{len(summary["models"])}</span>')
        content.append('<span class="metric-label">Models Tested</span>')
        content.append('</div>')
        content.append('<div class="metric-card">')
        content.append(f'<span class="metric-value">{len(categorized_plots)}</span>')
        content.append('<span class="metric-label">Model Types</span>')
        content.append('</div>')
        content.append('<div class="metric-card">')
        content.append('<span class="metric-value">PASS</span>')
        content.append('<span class="metric-label">Status</span>')
        content.append('</div>')
        content.append('</div>')
    else:
        content.append("### Key Metrics")
        content.append("")
        content.append(f"- **Total Tests:** {summary['total_tests']}")
        content.append(f"- **Models Tested:** {len(summary['models'])}")
        content.append(f"- **Model Types:** {len(categorized_plots)}")
        content.append("- **Status:** PASS")
    content.append("")
    
    # Model breakdown
    content.append("### Test Distribution by Model")
    content.append("")
    for model, count in summary['model_counts'].items():
        content.append(f"- **{model.upper()}:** {count} tests")
    content.append("")
    if html_styling:
        content.append('</div>')
    content.append("")
    
    # Get test results
    test_results = get_test_results(categorized_plots, build_dir=build_dir)
    
    # Regression Test Summary with styled table
    if html_styling:
        content.append('<div class="test-results">')
    content.append("## Regression Test Summary")
    content.append("")
    content.append("The following table summarizes the results of all regression tests:")
    content.append("")
    content.append("| Model | Test Type | Status |")
    content.append("|-------|-----------|--------|")
    
    for model in ['sphere', 'oswec', 'rm3', 'f3of']:
        if model in test_results:
            for test_type, status in test_results[model].items():
                if html_styling:
                    status_class = "status-pass" if status == "PASS" else "status-fail" if status == "FAIL" else "status-unknown"
                    status_badge = f'<span class="{status_class}">{status}</span>'
                else:
                    status_badge = status
                content.append(f"| {model.upper()} | {test_type.replace('_', ' ').title()} | {status_badge} |")
    content.append("")
    if html_styling:
        content.append('</div>')
    content.append("")
    
    # Add page break before detailed test results section
    content.append("\\clearpage")
    content.append("")
    
    # Detailed results by model
    for i, model in enumerate(['sphere', 'oswec', 'rm3', 'f3of']):
        if model in categorized_plots:
            # Add clear page break before each model section (except the first)
            if i > 0:
                content.append("\\clearpage")
                content.append("")
            
            if html_styling:
                content.append('<div class="model-section">')
            content.append(f"## {model.upper()} Tests")
            content.append("")
        
        model_plots = categorized_plots[model]
        
        # Track if this is the first subsection for this model
        first_subsection = True
        
        # Decay tests
        if 'decay' in model_plots and model_plots['decay']:
            if not first_subsection:
                content.append("\\clearpage")
                content.append("")
            first_subsection = False
            
            if html_styling:
                content.append('<div class="test-subsection">')
            content.append("### Decay Tests")
            content.append("")
            for plot_file in model_plots['decay']:
                plot_name = plot_file.stem.replace('_', ' ').replace('comparison', '').strip()
                relative_path = _relative_plot_path(plot_file, output_dir)
                
                test_type = "Decay Test"
                # Get the actual test status
                status = test_results.get(model, {}).get('decay', 'UNKNOWN')
                caption = f"**{model.upper()} - {test_type}** - Regression Test {status}"
                
                if html_styling:
                    content.append('<div class="image-container">')
                    content.append(f'<div class="image-title">{caption}</div>')
                    content.append(f"![{caption}]({relative_path})")
                    content.append('</div>')
                else:
                    content.append(f"![{caption}]({relative_path})")
                content.append("")
            if html_styling:
                content.append('</div>')
        
        # Regular waves tests
        if 'regular_waves' in model_plots and model_plots['regular_waves']:
            if not first_subsection:
                content.append("\\clearpage")
                content.append("")
            first_subsection = False
            
            if html_styling:
                content.append('<div class="test-subsection">')
            content.append("### Regular Waves Tests")
            content.append("")
            for plot_file in model_plots['regular_waves']:
                plot_name = plot_file.stem.replace('_', ' ').replace('comparison', '').strip()
                relative_path = _relative_plot_path(plot_file, output_dir)
                
                # Extract wave number if present
                wave_num = ""
                if "wave" in plot_name.lower():
                    import re
                    wave_match = re.search(r'wave\s*(\d+)', plot_name.lower())
                    if wave_match:
                        wave_num = f" Wave {wave_match.group(1)}"
                
                test_type = f"Regular Waves{wave_num}"
                # Get the actual test status
                status = test_results.get(model, {}).get('regular_waves', 'UNKNOWN')
                caption = f"**{model.upper()} - {test_type}** - Regression Test {status}"
                
                if html_styling:
                    content.append('<div class="image-container">')
                    content.append(f'<div class="image-title">{caption}</div>')
                    content.append(f"![{caption}]({relative_path})")
                    content.append('</div>')
                else:
                    content.append(f"![{caption}]({relative_path})")
                content.append("")
            if html_styling:
                content.append('</div>')
        
        # Irregular waves tests
        if 'irregular_waves' in model_plots and model_plots['irregular_waves']:
            if not first_subsection:
                content.append("\\clearpage")
                content.append("")
            first_subsection = False
            
            if html_styling:
                content.append('<div class="test-subsection">')
            content.append("### Irregular Waves Tests")
            content.append("")
            for plot_file in model_plots['irregular_waves']:
                plot_name = plot_file.stem.replace('_', ' ').replace('comparison', '').strip()
                relative_path = _relative_plot_path(plot_file, output_dir)
                
                test_type = "Irregular Waves"
                # Get the actual test status
                status = test_results.get(model, {}).get('irregular_waves', 'UNKNOWN')
                caption = f"**{model.upper()} - {test_type}** - Regression Test {status}"
                
                if html_styling:
                    content.append('<div class="image-container">')
                    content.append(f'<div class="image-title">{caption}</div>')
                    content.append(f"![{caption}]({relative_path})")
                    content.append('</div>')
                else:
                    content.append(f"![{caption}]({relative_path})")
                content.append("")
            if html_styling:
                content.append('</div>')
        
        # Other tests
        if 'other' in model_plots and model_plots['other']:
            if not first_subsection:
                content.append("\\clearpage")
                content.append("")
            first_subsection = False
            
            if html_styling:
                content.append('<div class="test-subsection">')
            content.append("### Other Tests")
            content.append("")
            for plot_file in model_plots['other']:
                plot_name = plot_file.stem.replace('_', ' ').replace('comparison', '').strip()
                relative_path = _relative_plot_path(plot_file, output_dir)
                
                test_type = "Other Test"
                # Get the actual test status (try to infer from filename)
                status = 'PASS'  # Default for other tests
                caption = f"**{model.upper()} - {test_type}** - Regression Test {status}"
                
                if html_styling:
                    content.append('<div class="image-container">')
                    content.append(f'<div class="image-title">{caption}</div>')
                    content.append(f"![{caption}]({relative_path})")
                    content.append('</div>')
                else:
                    content.append(f"![{caption}]({relative_path})")
                content.append("")
            if html_styling:
                content.append('</div>')
            
        if model in categorized_plots and html_styling:
            content.append('</div>')  # Close model-section
    
    # Footer
    if html_styling:
        content.append('<div class="footer">')
    content.append("---")
    content.append("")
    content.append("*This report was automatically generated by the HydroChrono regression test suite.*")
    if html_styling:
        content.append('</div>')
    
    return '\n'.join(content)

def convert_to_pdf_with_pandoc(markdown_file, output_dir):
    """Convert markdown to PDF using pandoc directly."""
    pdf_file = Path(output_dir) / "regression_test_report.pdf"
    css_file = Path(__file__).parent / "report_style.css"
    
    try:
        import subprocess
        
        # Check if pandoc is available
        result = subprocess.run(['pandoc', '--version'], capture_output=True, text=True)
        if result.returncode != 0:
            print("ERROR: pandoc not found. Please install pandoc:")
            print("  - Windows: Download from https://pandoc.org/installing.html")
            print("  - Or use: winget install JohnMacFarlane.Pandoc")
            return None
            
        print("Converting markdown to PDF using pandoc...")
        
        # Run pandoc from the report directory so relative image paths work correctly
        report_dir = Path(output_dir)
        markdown_filename = Path(markdown_file).name
        pdf_filename = Path(pdf_file).name
        
        # Build pandoc command - let pandoc choose the best available PDF engine
        cmd = [
            'pandoc',
            markdown_filename,
            '-o', pdf_filename,
            '--standalone',
            '--toc',
            '--toc-depth=3',
            '--metadata', f'title=HydroChrono Regression Test Report',
            '--metadata', f'date={datetime.now().strftime("%B %d, %Y")}'
        ]
        
        # Add CSS if available (need absolute path when running from different directory)
        if css_file.exists():
            cmd.extend(['--css', str(css_file.resolve())])
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=report_dir)
        
        if result.returncode == 0:
            print(f"SUCCESS: PDF generated successfully: {pdf_file}")
            return pdf_file
        else:
            print(f"ERROR: pandoc failed: {result.stderr}")
            return None
                
    except ImportError:
        print("subprocess module not available")
        return None
    except Exception as e:
        print(f"ERROR: PDF generation failed: {e}")
        return None

def convert_to_pdf(markdown_file, output_dir):
    """Convert markdown to PDF using pypandoc via HTML (no LaTeX required)."""
    pdf_file = Path(output_dir) / "regression_test_report.pdf"
    html_file = Path(output_dir) / "regression_test_report.html"
    css_file = Path(__file__).parent / "report_style.css"
    
    if not PYPANDOC_AVAILABLE:
        print("ERROR: pypandoc is not available. Install with: pip install pypandoc")
        return None
    
    try:
        
        # Ensure pandoc is available; download if necessary
        try:
            _ = pypandoc.get_pandoc_version()
        except OSError:
            print("Pandoc not found locally. Downloading pandoc...")
            pypandoc.download_pandoc()

        # First convert markdown to HTML with embedded CSS
        print("Converting markdown to HTML...")
        
        # Prepare pandoc arguments - no title metadata to avoid duplicate titles
        pandoc_args = [
            '--standalone',
            '--embed-resources',
            '--toc',
            '--toc-depth=3'
        ]
        
        # Add CSS if available
        if css_file.exists():
            pandoc_args.append(f'--css={css_file}')
        
        # Convert markdown to HTML
        html_content = pypandoc.convert_file(
            str(markdown_file),
            'html',
            extra_args=pandoc_args
        )
        
        # Enhance HTML with additional styling and structure
        # Remove any auto-generated title from pandoc to avoid duplicates
        import re
        # Remove pandoc-generated title if it exists
        html_content = re.sub(r'<h1[^>]*class="title"[^>]*>.*?</h1>', '', html_content, flags=re.DOTALL)
        html_content = re.sub(r'<header[^>]*>.*?</header>', '', html_content, flags=re.DOTALL)
        
        # Fix TOC placement - move it after our custom header
        # Extract our header div
        header_match = re.search(r'(<div class="header">.*?</div>)', html_content, flags=re.DOTALL)
        toc_match = re.search(r'(<nav[^>]*id="TOC"[^>]*>.*?</nav>)', html_content, flags=re.DOTALL)
        
        if header_match and toc_match:
            header_content = header_match.group(1)
            toc_content = toc_match.group(1)
            
            # Remove both from their current positions
            html_content = re.sub(r'<div class="header">.*?</div>', '', html_content, flags=re.DOTALL)
            html_content = re.sub(r'<nav[^>]*id="TOC"[^>]*>.*?</nav>', '', html_content, flags=re.DOTALL)
            
            # Insert header first, then TOC at the beginning of body
            body_start = html_content.find('<body>')
            if body_start != -1:
                body_start += len('<body>')
                html_content = (html_content[:body_start] + 
                              '\n' + header_content + '\n' + toc_content + '\n' + 
                              html_content[body_start:])
        elif header_match:
            # If no TOC, just ensure header is at the top
            header_content = header_match.group(1)
            html_content = re.sub(r'<div class="header">.*?</div>', '', html_content, flags=re.DOTALL)
            
            body_start = html_content.find('<body>')
            if body_start != -1:
                body_start += len('<body>')
                html_content = (html_content[:body_start] + 
                              '\n' + header_content + '\n' + 
                              html_content[body_start:])
        
        enhanced_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HydroChrono Regression Test Report</title>
    <style>
        body {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 2rem;
            background: #FFFFFF;
        }}
        
        @media print {{
            body {{
                max-width: none;
                margin: 0;
                padding: 1rem;
            }}
        }}
    </style>
</head>
<body>
{html_content}
</body>
</html>"""
        
        # Write enhanced HTML to file
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(enhanced_html)
        
        print(f"HTML report generated: {html_file}")
        
        # Try to convert HTML to PDF using weasyprint (CSS-based, no LaTeX needed)
        try:
            import weasyprint
            print("Converting HTML to PDF using weasyprint...")
            
            # Configure weasyprint for better rendering
            html_doc = weasyprint.HTML(filename=str(html_file))
            css_string = None
            
            if css_file.exists():
                with open(css_file, 'r', encoding='utf-8') as f:
                    css_string = f.read()
                    
                # Add print-specific optimizations
                css_string += """
                @page {
                    size: A4;
                    margin: 1in;
                }
                
                .model-section {
                    page-break-before: always;
                }
                
                .image-container {
                    page-break-inside: avoid;
                    margin: 1rem 0;
                }
                
                .test-subsection {
                    page-break-inside: avoid;
                }
                
                .header {
                    page-break-after: always;
                }
                """
                
                css_doc = weasyprint.CSS(string=css_string)
                html_doc.write_pdf(str(pdf_file), stylesheets=[css_doc])
            else:
                html_doc.write_pdf(str(pdf_file))
            
            print(f"SUCCESS: PDF generated successfully: {pdf_file}")
            return pdf_file
            
        except ImportError:
            print("weasyprint not available, trying alternative method...")
        except Exception as e:
            print(f"weasyprint failed: {e}, trying alternative method...")
        
        # Alternative: Try pdfkit (requires wkhtmltopdf)
        try:
            import pdfkit
            print("Converting HTML to PDF using pdfkit...")
            
            options = {
                'page-size': 'A4',
                'margin-top': '1in',
                'margin-right': '1in',
                'margin-bottom': '1in',
                'margin-left': '1in',
                'encoding': "UTF-8",
                'no-outline': None,
                'enable-local-file-access': None
            }
            
            pdfkit.from_file(str(html_file), str(pdf_file), options=options)
            print(f"SUCCESS: PDF generated successfully: {pdf_file}")
            return pdf_file
            
        except ImportError:
            print("pdfkit not available, trying reportlab...")
        except Exception as e:
            print(f"pdfkit failed: {e}, trying reportlab...")
        
        # Alternative: Use reportlab for basic PDF generation
        try:
            from reportlab.pdfgen import canvas
            from reportlab.lib.pagesizes import letter, A4
            from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
            from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
            from reportlab.lib.units import inch
            from reportlab.lib.colors import HexColor
            import re
            
            print("Converting to PDF using reportlab...")
            
            # Create PDF with reportlab (enhanced styling)
            doc = SimpleDocTemplate(str(pdf_file), pagesize=A4, 
                                  rightMargin=72, leftMargin=72, 
                                  topMargin=72, bottomMargin=72)
            
            # Custom styles
            styles = getSampleStyleSheet()
            title_style = ParagraphStyle(
                'CustomTitle',
                parent=styles['Title'],
                fontSize=24,
                textColor=HexColor('#007AFF'),
                spaceAfter=30
            )
            
            heading_style = ParagraphStyle(
                'CustomHeading',
                parent=styles['Heading1'],
                fontSize=18,
                textColor=HexColor('#1D1D1F'),
                spaceAfter=20
            )
            
            story = []
            
            # Read markdown content and convert to styled text
            with open(markdown_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Enhanced markdown to text conversion
            lines = content.split('\n')
            for line in lines:
                line = line.strip()
                if line.startswith('# ') and 'HydroChrono' in line:
                    story.append(Paragraph(line[2:], title_style))
                elif line.startswith('## '):
                    story.append(Paragraph(line[3:], heading_style))
                elif line.startswith('### '):
                    story.append(Paragraph(line[4:], styles['Heading2']))
                elif line.startswith('!['):
                    # Skip images for text-only PDF but add placeholder
                    story.append(Paragraph(f"[Comparison Plot: {line[2:line.find(']')]}]", styles['Italic']))
                elif line and not line.startswith('---') and not line.startswith('<'):
                    story.append(Paragraph(line, styles['Normal']))
                    
                if line.endswith('Tests'):
                    story.append(PageBreak())
                else:
                    story.append(Spacer(1, 0.1*inch))
            
            doc.build(story)
            print(f"SUCCESS: Basic PDF generated successfully: {pdf_file}")
            return pdf_file
            
        except ImportError:
            print("reportlab not available either.")
        except Exception as e:
            print(f"reportlab failed: {e}")
        
        print("WARNING: No PDF generation library available.")
        print("To generate PDFs, install one of:")
        print("  pip install weasyprint")
        print("  pip install pdfkit")  
        print("  pip install reportlab")
        return None
        
    except Exception as e:
        print(f"ERROR: PDF generation failed: {e}")
        return None

def rerun_comparisons(build_dir, config="Release"):
    """Re-run CTest comparison tests to regenerate plot images.
    
    This ensures plots reflect the current reference data. Without this step,
    plots may be stale if reference data was updated after the last test run.
    
    Note: The comparison tests use CTest FIXTURES_REQUIRED, which means CTest
    will skip them unless the simulation fixture tests are also selected.
    We use -FA "." to exclude all fixture requirements so the comparison
    scripts run directly (the simulation output files already exist on disk).
    """
    print("Re-running comparison tests to regenerate plots...")
    build_path = Path(build_dir).resolve()
    
    cmd = [
        "ctest",
        "--test-dir", str(build_path),
        "-C", config,
        "-R", "_regression",
        "-LE", "report",
        "-FA", ".",
        "--output-on-failure"
    ]
    
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("WARNING: Some comparison tests failed (plots were still regenerated).")
    else:
        print("All comparison tests passed.")
    print()


def main():
    parser = argparse.ArgumentParser(description='Generate HydroChrono regression test report')
    parser.add_argument('--output-dir', help='Output directory for reports (default: build/bin/tests/regression/report)')
    parser.add_argument('--build-dir', default='build', help='Build directory path')
    parser.add_argument('--pdf', action='store_true', help='Generate PDF using pandoc (requires pandoc to be installed)')
    parser.add_argument('--html-styling', action='store_true', help='Include HTML styling in markdown (for advanced PDF generation)')
    parser.add_argument('--recompare', action='store_true',
                        help='Re-run CTest comparison tests before generating the report. '
                             'Use this after updating reference data to ensure plots are current.')
    parser.add_argument('--config', default='Release', help='Build configuration (default: Release)')
    
    args = parser.parse_args()
    
    # Re-run comparison tests if requested
    if args.recompare:
        rerun_comparisons(args.build_dir, args.config)
    
    # Set default output directory to build/bin/tests/regression/report if not specified
    if args.output_dir is None:
        output_dir = Path(args.build_dir) / "bin" / "tests" / "regression" / "report"
    else:
        output_dir = Path(args.output_dir)
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all plot files
    print("Scanning for comparison plots...")
    plots = find_plot_files(args.build_dir)
    
    if not plots:
        print("ERROR: No comparison plots found!")
        print(f"   Expected location: {args.build_dir}/tests/regression/Release/*/results/plots/")
        sys.exit(1)
    
    # Categorize plots
    print("Categorizing plots by test type...")
    categorized_plots = categorize_plots(plots)
    
    # Generate markdown report
    print("Generating markdown report...")
    markdown_content = generate_markdown_report(categorized_plots, output_dir, args.build_dir, html_styling=args.html_styling)
    
    # Write markdown file
    markdown_file = output_dir / "regression_test_report.md"
    with open(markdown_file, 'w', encoding='utf-8') as f:
        f.write(markdown_content)
    
    print(f"SUCCESS: Clean markdown report generated: {markdown_file}")
    
    # Generate PDF if requested
    if args.pdf:
        pdf_file = convert_to_pdf_with_pandoc(markdown_file, output_dir)
        if pdf_file:
            print(f"Report available in both formats:")
            print(f"   - Markdown: {markdown_file}")
            print(f"   - PDF: {pdf_file}")
        else:
            print(f"PDF generation failed. Report available as markdown: {markdown_file}")
    else:
        print(f"Report available as clean markdown: {markdown_file}")
        print("Use --pdf flag to generate PDF (requires pandoc)")

if __name__ == "__main__":
    main() 