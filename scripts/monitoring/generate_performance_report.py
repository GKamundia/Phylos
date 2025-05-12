#!/usr/bin/env python3
"""
Generate performance report from benchmark files
"""

import os
import sys
import json
import time
from datetime import datetime
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import snakemake

# Add parent directory to path for importing utils
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.log_utils import setup_logger, log_execution_stats, log_with_context

# Configure logger
logger = setup_logger(
    name="generate_performance_report", 
    log_file=snakemake.log[0] if "snakemake" in globals() else "generate_performance_report.log", 
    level="INFO"
)

def parse_benchmark_file(file_path):
    """Parse Snakemake benchmark file into a DataFrame"""
    try:
        # Load benchmark data
        df = pd.read_csv(file_path, sep='\t')
        
        # Extract rule name from filename
        rule_name = file_path.stem
        # Extract rule name without prefix/suffix
        if '_' in rule_name:
            parts = rule_name.split('_')
            if len(parts) >= 2:
                rule_name = '_'.join(parts[:-1])  # Remove last part (often a timestamp or suffix)
        
        # Add rule name to DataFrame
        df['rule'] = rule_name
        
        return df
    except Exception as e:
        logger.error(f"Error parsing benchmark file {file_path}: {e}")
        return None

def generate_json_report(benchmarks, output_file):
    """Generate JSON performance report"""
    # Combine all benchmark data
    all_data = []
    for rule, df in benchmarks.items():
        if not df.empty:
            stats = {
                "rule": rule,
                "runtime_s": float(df['s'].mean()),
                "max_runtime_s": float(df['s'].max()),
                "min_runtime_s": float(df['s'].min()),
                "mean_mem_mb": float(df['max_rss'] / 1024).mean() if 'max_rss' in df else None,
                "max_mem_mb": float(df['max_rss'] / 1024).max() if 'max_rss' in df else None,
            }
            all_data.append(stats)
    
    # Sort by runtime (descending)
    all_data.sort(key=lambda x: x['runtime_s'], reverse=True)
    
    # Create report
    report = {
        "generated_at": datetime.now().isoformat(),
        "pathogen": snakemake.params.pathogen,
        "pathogen_name": snakemake.params.pathogen_name,
        "segment_mode": snakemake.params.segment_mode,
        "segments": snakemake.params.segments,
        "metrics": {
            "rules": all_data,
            "total_runtime_s": sum(item['runtime_s'] for item in all_data),
            "max_runtime_rule": max(all_data, key=lambda x: x['runtime_s'])['rule'] if all_data else None,
            "max_memory_rule": max(all_data, key=lambda x: x['max_mem_mb'] if x['max_mem_mb'] is not None else 0)['rule'] if all_data else None
        }
    }
    
    # Write report
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    return report

def generate_html_dashboard(benchmarks, output_file, report):
    """Generate HTML performance dashboard with visualizations"""
    # Create a more detailed DataFrame for visualization
    data = []
    for rule, df in benchmarks.items():
        if not df.empty:
            for _, row in df.iterrows():
                entry = {
                    'rule': rule,
                    'runtime_s': row['s'],
                    'max_rss_mb': row['max_rss'] / 1024 if 'max_rss' in df else 0,
                    'io_in_mb': row['io_in'] / (1024 * 1024) if 'io_in' in df else 0,
                    'io_out_mb': row['io_out'] / (1024 * 1024) if 'io_out' in df else 0
                }
                data.append(entry)
    
    if not data:
        logger.warning("No benchmark data available for visualization")
        with open(output_file, 'w') as f:
            f.write("<html><body><h1>No benchmark data available</h1></body></html>")
        return
    
    df = pd.DataFrame(data)
    
    # Create visualization
    plt.figure(figsize=(12, 10))
    
    # Runtime by rule
    plt.subplot(2, 1, 1)
    runtime_by_rule = df.groupby('rule')['runtime_s'].mean().sort_values(ascending=False)
    sns.barplot(x=runtime_by_rule.values, y=runtime_by_rule.index)
    plt.title(f'Average Runtime by Rule ({snakemake.params.pathogen_name})')
    plt.xlabel('Runtime (seconds)')
    plt.tight_layout()
    
    # Memory usage by rule
    plt.subplot(2, 1, 2)
    if 'max_rss_mb' in df and df['max_rss_mb'].sum() > 0:
        mem_by_rule = df.groupby('rule')['max_rss_mb'].mean().sort_values(ascending=False)
        sns.barplot(x=mem_by_rule.values, y=mem_by_rule.index)
        plt.title('Average Memory Usage by Rule')
        plt.xlabel('Memory (MB)')
    else:
        plt.text(0.5, 0.5, "No memory usage data available", ha='center')
    
    plt.tight_layout()
    
    # Save visualization
    plt.savefig('temp_performance_viz.png')
    
    # Create HTML
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Nextstrain Performance Dashboard - {snakemake.params.pathogen_name}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            .container {{ max-width: 1200px; margin: 0 auto; }}
            .header {{ background-color: #f8f9fa; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
            .chart {{ margin-top: 30px; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ text-align: left; padding: 8px; border-bottom: 1px solid #ddd; }}
            th {{ background-color: #f2f2f2; }}
            tr:hover {{ background-color: #f5f5f5; }}
            .highlight {{ background-color: #ffffcc; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>Nextstrain Performance Dashboard</h1>
                <p>Pathogen: <strong>{snakemake.params.pathogen_name}</strong> ({snakemake.params.pathogen})</p>
                <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p>Segment Mode: {snakemake.params.segment_mode}</p>
            </div>
            
            <h2>Performance Summary</h2>
            <p><strong>Total Runtime:</strong> {report['metrics']['total_runtime_s']:.2f} seconds ({report['metrics']['total_runtime_s']/60:.2f} minutes)</p>
            <p><strong>Longest Running Rule:</strong> {report['metrics']['max_runtime_rule'] or 'N/A'}</p>
            <p><strong>Most Memory-Intensive Rule:</strong> {report['metrics']['max_memory_rule'] or 'N/A'}</p>
            
            <div class="chart">
                <img src="data:image/png;base64,{open('temp_performance_viz.png', 'rb').read().hex()}" width="100%">
            </div>
            
            <h2>Detailed Rule Performance</h2>
            <table>
                <tr>
                    <th>Rule</th>
                    <th>Avg Runtime (s)</th>
                    <th>Max Runtime (s)</th>
                    <th>Min Runtime (s)</th>
                    <th>Avg Memory (MB)</th>
                    <th>Max Memory (MB)</th>
                </tr>
    """
    
    # Add rows for each rule
    for rule_data in sorted(report['metrics']['rules'], key=lambda x: x['runtime_s'], reverse=True):
        highlight = ' class="highlight"' if rule_data['rule'] in [report['metrics']['max_runtime_rule'], report['metrics']['max_memory_rule']] else ''
        html += f"""
                <tr{highlight}>
                    <td>{rule_data['rule']}</td>
                    <td>{rule_data['runtime_s']:.2f}</td>
                    <td>{rule_data['max_runtime_s']:.2f}</td>
                    <td>{rule_data['min_runtime_s']:.2f}</td>
                    <td>{rule_data['mean_mem_mb']:.2f if rule_data['mean_mem_mb'] is not None else 'N/A'}</td>
                    <td>{rule_data['max_mem_mb']:.2f if rule_data['max_mem_mb'] is not None else 'N/A'}</td>
                </tr>
        """
    
    html += """
            </table>
        </div>
    </body>
    </html>
    """
    
    # Write HTML file
    with open(output_file, 'w') as f:
        f.write(html)
    
    # Clean up
    if os.path.exists('temp_performance_viz.png'):
        os.remove('temp_performance_viz.png')

def main():
    start_time = time.time()
    
    try:
        # Get benchmark files
        benchmark_files = snakemake.input.benchmarks
        logger.info(f"Processing {len(benchmark_files)} benchmark files")
        
        # Parse benchmark files
        benchmarks = {}
        for file_path in benchmark_files:
            file_path = Path(file_path)
            df = parse_benchmark_file(file_path)
            if df is not None:
                rule_name = df['rule'].iloc[0]
                if rule_name in benchmarks:
                    benchmarks[rule_name] = pd.concat([benchmarks[rule_name], df])
                else:
                    benchmarks[rule_name] = df
        
        logger.info(f"Successfully parsed {len(benchmarks)} rules' benchmark data")
        
        # Generate JSON report
        report = generate_json_report(benchmarks, snakemake.output.report)
        logger.info(f"Generated JSON performance report at {snakemake.output.report}")
        
        # Generate HTML dashboard
        generate_html_dashboard(benchmarks, snakemake.output.html, report)
        logger.info(f"Generated HTML performance dashboard at {snakemake.output.html}")
        
        # Log execution statistics
        log_execution_stats(logger, start_time, {
            "operation": "generate_performance_report",
            "rule_count": len(benchmarks),
            "output_files": [str(snakemake.output.report), str(snakemake.output.html)]
        })
    
    except Exception as e:
        logger.error(f"Error generating performance report: {str(e)}")
        log_execution_stats(logger, start_time, {"operation": "generate_performance_report", "error": str(e)}, status="failed")
        raise

if __name__ == "__main__":
    if "snakemake" not in globals():
        logger.error("This script should be run via Snakemake")
        sys.exit(1)
    main()