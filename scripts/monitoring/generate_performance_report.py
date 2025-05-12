#!/usr/bin/env python3
"""
Generate performance reports and visualizations from benchmarking data
"""

import os
import sys
import json
import time
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional

import snakemake

# Add parent directory to path for importing utils
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.log_utils import setup_logger, log_with_context, log_execution_stats

# Configure logger
logger = setup_logger(
    name="performance_report",
    log_file=snakemake.log[0] if "snakemake" in globals() else "generate_performance_report.log"
)

def parse_benchmark_files(benchmark_files: List[str]) -> pd.DataFrame:
    """Parse Snakemake benchmark files into a DataFrame"""
    all_data = []
    
    for file_path in benchmark_files:
        try:
            # Extract rule name from the benchmark filename
            file_name = os.path.basename(file_path)
            rule_name = file_name.split("_")[0]
            
            # Read benchmark data
            df = pd.read_csv(file_path, sep="\t")
            
            # Add rule name and file path
            df["rule"] = rule_name
            df["benchmark_file"] = file_path
            
            all_data.append(df)
        except Exception as e:
            logger.error(f"Error parsing benchmark file {file_path}: {e}")
    
    # Combine all data
    if all_data:
        result = pd.concat(all_data, ignore_index=True)
        return result
    else:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=["s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss", "io_in", "io_out", "mean_load", "cpu_time", "rule", "benchmark_file"])

def generate_summary(benchmark_data: pd.DataFrame) -> Dict[str, Any]:
    """Generate summary statistics from benchmark data"""
    if benchmark_data.empty:
        return {
            "total_rules": 0,
            "total_runtime_seconds": 0,
            "total_runtime_formatted": "0:00:00",
            "most_time_consuming_rule": None,
            "highest_memory_rule": None,
            "rules_summary": []
        }
    
    # Calculate summary statistics
    rules_summary = []
    
    # Group by rule name
    grouped = benchmark_data.groupby("rule")
    
    for rule, group in grouped:
        # Calculate statistics for this rule
        runtime_seconds = group["s"].sum()
        max_memory_mb = group["max_rss"].max() / 1024  # Convert to MB
        
        # Format runtime
        hours, remainder = divmod(runtime_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        runtime_formatted = f"{int(hours)}:{int(minutes):02d}:{seconds:.2f}"
        
        rules_summary.append({
            "rule": rule,
            "runtime_seconds": runtime_seconds,
            "runtime_formatted": runtime_formatted,
            "max_memory_mb": max_memory_mb,
            "cpu_time": group["cpu_time"].sum(),
            "instances": len(group)
        })
    
    # Sort rules by runtime
    rules_summary.sort(key=lambda x: x["runtime_seconds"], reverse=True)
    
    # Calculate total runtime
    total_runtime_seconds = benchmark_data["s"].sum()
    hours, remainder = divmod(total_runtime_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    total_runtime_formatted = f"{int(hours)}:{int(minutes):02d}:{seconds:.2f}"
    
    # Find most time-consuming rule
    most_time_consuming_rule = rules_summary[0] if rules_summary else None
    
    # Find highest memory rule
    highest_memory_rule = max(rules_summary, key=lambda x: x["max_memory_mb"]) if rules_summary else None
    
    return {
        "total_rules": len(rules_summary),
        "total_runtime_seconds": total_runtime_seconds,
        "total_runtime_formatted": total_runtime_formatted,
        "most_time_consuming_rule": most_time_consuming_rule,
        "highest_memory_rule": highest_memory_rule,
        "rules_summary": rules_summary
    }

def generate_performance_plots(benchmark_data: pd.DataFrame, output_dir: str, prefix: str) -> List[str]:
    """Generate performance visualization plots"""
    if benchmark_data.empty:
        return []
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    generated_files = []
    
    try:
        # 1. Runtime by rule
        plt.figure(figsize=(10, 8))
        rule_times = benchmark_data.groupby("rule")["s"].sum().sort_values(ascending=False)
        rule_times.plot(kind="bar")
        plt.title("Runtime by Rule")
        plt.ylabel("Time (seconds)")
        plt.xlabel("Rule")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        runtime_plot = os.path.join(output_dir, f"{prefix}_runtime_by_rule.png")
        plt.savefig(runtime_plot)
        plt.close()
        generated_files.append(runtime_plot)
        
        # 2. Memory usage by rule
        plt.figure(figsize=(10, 8))
        rule_memory = benchmark_data.groupby("rule")["max_rss"].max().sort_values(ascending=False) / 1024  # Convert to MB
        rule_memory.plot(kind="bar")
        plt.title("Maximum Memory Usage by Rule")
        plt.ylabel("Memory (MB)")
        plt.xlabel("Rule")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        memory_plot = os.path.join(output_dir, f"{prefix}_memory_by_rule.png")
        plt.savefig(memory_plot)
        plt.close()
        generated_files.append(memory_plot)
        
        # 3. CPU time vs. wall time
        plt.figure(figsize=(10, 8))
        benchmark_data["cpu_efficiency"] = (benchmark_data["cpu_time"] / benchmark_data["s"]) * 100
        rule_efficiency = benchmark_data.groupby("rule")["cpu_efficiency"].mean().sort_values(ascending=False)
        rule_efficiency.plot(kind="bar")
        plt.title("CPU Efficiency by Rule (CPU Time / Wall Time)")
        plt.ylabel("Efficiency (%)")
        plt.xlabel("Rule")
        plt.xticks(rotation=45, ha="right")
        plt.axhline(y=100, color="r", linestyle="--", label="Ideal (100%)")
        plt.legend()
        plt.tight_layout()
        efficiency_plot = os.path.join(output_dir, f"{prefix}_cpu_efficiency.png")
        plt.savefig(efficiency_plot)
        plt.close()
        generated_files.append(efficiency_plot)
        
    except Exception as e:
        logger.error(f"Error generating performance plots: {e}")
    
    return generated_files

def generate_html_report(summary: Dict[str, Any], plots: List[str], output_file: str, params: Dict[str, Any]) -> bool:
    """Generate an HTML dashboard with performance data"""
    try:
        # Convert plot paths to relative paths for HTML
        plot_basenames = [os.path.basename(plot) for plot in plots]
        
        html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Nextstrain Performance Report</title>
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    line-height: 1.6;
                    margin: 0;
                    padding: 20px;
                    color: #333;
                }}
                h1, h2, h3 {{
                    color: #2c3e50;
                }}
                .container {{
                    max-width: 1200px;
                    margin: 0 auto;
                }}
                .summary-card {{
                    background-color: #f8f9fa;
                    border-radius: 5px;
                    padding: 20px;
                    margin-bottom: 20px;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }}
                .metric {{
                    display: inline-block;
                    margin-right: 30px;
                    margin-bottom: 15px;
                }}
                .metric-value {{
                    font-size: 24px;
                    font-weight: bold;
                    display: block;
                }}
                .metric-label {{
                    font-size: 14px;
                    color: #666;
                }}
                table {{
                    width: 100%;
                    border-collapse: collapse;
                    margin: 20px 0;
                }}
                th, td {{
                    border: 1px solid #ddd;
                    padding: 8px;
                    text-align: left;
                }}
                th {{
                    background-color: #f2f2f2;
                }}
                tr:nth-child(even) {{
                    background-color: #f9f9f9;
                }}
                .plot-container {{
                    margin: 20px 0;
                    text-align: center;
                }}
                .plot-container img {{
                    max-width: 100%;
                    height: auto;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>Nextstrain Performance Report</h1>
                <p>
                    <strong>Pathogen:</strong> {params.get('pathogen_name', 'Unknown')} ({params.get('pathogen', 'unknown')})<br>
                    <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br>
                    <strong>Segment Mode:</strong> {params.get('segment_mode', 'unknown')}
                </p>
                
                <div class="summary-card">
                    <h2>Performance Summary</h2>
                    <div class="metric">
                        <span class="metric-value">{summary['total_rules']}</span>
                        <span class="metric-label">Total Rules</span>
                    </div>
                    <div class="metric">
                        <span class="metric-value">{summary['total_runtime_formatted']}</span>
                        <span class="metric-label">Total Runtime</span>
                    </div>
        """
        
        # Add most time consuming rule if available
        if summary.get('most_time_consuming_rule'):
            html_content += f"""
                    <div class="metric">
                        <span class="metric-value">{summary['most_time_consuming_rule']['rule']}</span>
                        <span class="metric-label">Most Time-Consuming Rule</span>
                    </div>
                    <div class="metric">
                        <span class="metric-value">{summary['most_time_consuming_rule']['runtime_formatted']}</span>
                        <span class="metric-label">Runtime of Slowest Rule</span>
                    </div>
            """
        
        # Add highest memory rule if available
        if summary.get('highest_memory_rule'):
            html_content += f"""
                    <div class="metric">
                        <span class="metric-value">{summary['highest_memory_rule']['rule']}</span>
                        <span class="metric-label">Highest Memory Rule</span>
                    </div>
                    <div class="metric">
                        <span class="metric-value">{summary['highest_memory_rule']['max_memory_mb']:.2f} MB</span>
                        <span class="metric-label">Max Memory</span>
                    </div>
            """
        
        html_content += """
                </div>
                
                <h2>Rule Performance</h2>
                <table>
                    <tr>
                        <th>Rule</th>
                        <th>Runtime</th>
                        <th>Max Memory (MB)</th>
                        <th>CPU Time (s)</th>
                        <th>Instances</th>
                    </tr>
        """
        
        # Add rows for each rule
        for rule in summary.get('rules_summary', []):
            html_content += f"""
                    <tr>
                        <td>{rule['rule']}</td>
                        <td>{rule['runtime_formatted']}</td>
                        <td>{rule['max_memory_mb']:.2f}</td>
                        <td>{rule['cpu_time']:.2f}</td>
                        <td>{rule['instances']}</td>
                    </tr>
            """
        
        html_content += """
                </table>
                
                <h2>Performance Visualizations</h2>
        """
        
        # Add plots
        for plot in plot_basenames:
            html_content += f"""
                <div class="plot-container">
                    <img src="{plot}" alt="Performance Plot">
                </div>
            """
        
        html_content += """
            </div>
        </body>
        </html>
        """
        
        # Write HTML to file
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        return True
    except Exception as e:
        logger.error(f"Error generating HTML report: {e}")
        return False

def main():
    start_time = time.time()
    
    try:
        # Get parameters from Snakemake
        benchmark_files = snakemake.input.benchmarks
        output_report_json = snakemake.output.report
        output_report_html = snakemake.output.html
        params = {
            "pathogen": snakemake.params.pathogen,
            "pathogen_name": snakemake.params.pathogen_name,
            "segment_mode": snakemake.params.segment_mode,
            "segments": snakemake.params.segments
        }
        
        # Create output directories
        os.makedirs(os.path.dirname(output_report_json), exist_ok=True)
        os.makedirs(os.path.dirname(output_report_html), exist_ok=True)
        
        # Convert to strings if Path objects
        benchmark_files = [str(f) for f in benchmark_files]
        
        # Parse benchmark data
        log_with_context(logger, "INFO", f"Parsing {len(benchmark_files)} benchmark files")
        benchmark_data = parse_benchmark_files(benchmark_files)
        
        # Generate summary
        log_with_context(logger, "INFO", "Generating performance summary")
        summary = generate_summary(benchmark_data)
        
        # Generate plots
        log_with_context(logger, "INFO", "Generating performance plots")
        plots_dir = os.path.dirname(output_report_html)
        prefix = os.path.splitext(os.path.basename(output_report_html))[0]
        plots = generate_performance_plots(benchmark_data, plots_dir, prefix)
        
        # Generate HTML report
        log_with_context(logger, "INFO", "Generating HTML report")
        html_generated = generate_html_report(summary, plots, output_report_html, params)
        
        # Write JSON report
        with open(output_report_json, 'w') as f:
            json.dump(summary, f, indent=2)
        
        log_execution_stats(logger, start_time, {
            "operation": "generate_performance_report",
            "benchmark_files": len(benchmark_files),
            "rules_analyzed": summary["total_rules"],
            "plots_generated": len(plots),
            "html_generated": html_generated
        })
        
        return 0
    
    except Exception as e:
        log_with_context(logger, "ERROR", f"Error generating performance report: {e}")
        log_execution_stats(logger, start_time, {
            "operation": "generate_performance_report",
            "error": str(e)
        }, status="failed")
        return 1

if __name__ == "__main__":
    if "snakemake" not in globals():
        logger.error("This script is designed to be run through Snakemake")
        sys.exit(1)
    
    sys.exit(main())