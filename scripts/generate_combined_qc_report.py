#!/usr/bin/env python3
"""
Generate a combined QC report for all segments
"""

import json
import logging

import snakemake

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get input parameters from Snakemake
segment_qc_files = snakemake.input.segment_qc
combined_report = snakemake.output.combined_report
pathogen = snakemake.params.pathogen
segments = snakemake.params.segments
log_file = snakemake.log[0]

# Configure logging to file
file_handler = logging.FileHandler(log_file)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(file_handler)

# Load QC data from each segment
segment_qc_data = {}
for segment, qc_file in zip(segments, segment_qc_files):
    with open(qc_file, 'r') as f:
        segment_qc_data[segment] = json.load(f)
    logger.info(f"Loaded QC data for segment {segment}")

# Generate combined HTML report
html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>{pathogen} Multi-Segment QC Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1, h2, h3 {{ color: #2c3e50; }}
        .summary {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        .segment {{ margin-bottom: 30px; border-left: 5px solid #2c3e50; padding-left: 15px; }}
        .metric {{ margin-bottom: 10px; }}
        .metric-name {{ font-weight: bold; }}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        .good {{ color: green; }}
        .warning {{ color: orange; }}
        .error {{ color: red; }}
    </style>
</head>
<body>
    <h1>{pathogen} Multi-Segment QC Report</h1>
    
    <div class="summary">
        <h2>Summary</h2>
        <table>
            <tr>
                <th>Segment</th>
                <th>Initial Sequences</th>
                <th>Filtered Sequences</th>
                <th>QC Passed</th>
                <th>Pass Rate</th>
            </tr>
"""

# Add summary rows for each segment
for segment in segments:
    qc_data = segment_qc_data[segment]
    raw_count = qc_data.get("sequence_counts", {}).get("raw", 0)
    filtered_count = qc_data.get("sequence_counts", {}).get("filtered", 0)
    passed_count = qc_data.get("sequence_counts", {}).get("nextclade_passed", 0)
    pass_rate = f"{passed_count / raw_count * 100:.1f}%" if raw_count > 0 else "N/A"
    
    html_content += f"""
            <tr>
                <td>{segment}</td>
                <td>{raw_count}</td>
                <td>{filtered_count}</td>
                <td>{passed_count}</td>
                <td>{pass_rate}</td>
            </tr>
    """

# Close summary table and add segment-specific details
html_content += """
        </table>
    </div>
"""

# Add detailed section for each segment
for segment in segments:
    qc_data = segment_qc_data[segment]
    html_content += f"""
    <div class="segment">
        <h2>Segment {segment} Details</h2>
        
        <h3>Filtering Results</h3>
        <table>
            <tr>
                <th>Filter Type</th>
                <th>Sequences Excluded</th>
                <th>Threshold</th>
            </tr>
    """
    
    # Add filter results
    for filter_type, details in qc_data.get('filter_reasons', {}).items():
        if isinstance(details, dict) and 'excluded' in details:
            html_content += f"""
            <tr>
                <td>{filter_type.capitalize()}</td>
                <td>{details['excluded']}</td>
                <td>{details.get('threshold', 'N/A')}</td>
            </tr>
            """
    
    html_content += """
        </table>
        
        <h3>Nextclade QC Results</h3>
    """
    
    # Add Nextclade results if available
    if 'nextclade_stats' in qc_data and 'qc_scores' in qc_data['nextclade_stats']:
        html_content += """
        <table>
            <tr>
                <th>QC Score Range</th>
                <th>Count</th>
            </tr>
        """
        
        for score_range, count in sorted(qc_data['nextclade_stats']['qc_scores'].items()):
            html_content += f"""
            <tr>
                <td>{score_range}</td>
                <td>{count}</td>
            </tr>
            """
        
        html_content += """
        </table>
        """
    else:
        html_content += "<p>No Nextclade QC data available for this segment.</p>"
    
    html_content += """
    </div>
    """

# Close the HTML document
html_content += """
</body>
</html>
"""

# Write the HTML report
with open(combined_report, 'w') as f:
    f.write(html_content)

logger.info(f"Combined QC report written to {combined_report}")