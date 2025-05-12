#!/usr/bin/env python3
"""
Generates comprehensive QC reports for sequence data processing.
"""

import os
import sys
import json
import argparse
import logging
from datetime import datetime
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import re

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate comprehensive QC reports")
    
    # Input files
    parser.add_argument('--raw-sequences', required=True,
                        help="Path to raw input sequences FASTA")
    parser.add_argument('--filtered-sequences', required=True,
                        help="Path to filtered sequences FASTA")
    parser.add_argument('--nextclade-json', required=True,
                        help="Path to Nextclade JSON output")
    parser.add_argument('--nextclade-passed', required=True,
                        help="Path to Nextclade-filtered sequences FASTA")
    parser.add_argument('--raw-metadata', required=True,
                        help="Path to raw metadata TSV")
    parser.add_argument('--filtered-metadata', required=True,
                        help="Path to filtered metadata TSV")
    parser.add_argument('--filter-log', required=True,
                        help="Path to augur filter log file")
    parser.add_argument('--nextclade-log', required=True,
                        help="Path to Nextclade filtering log file")
    
    # Output files
    parser.add_argument('--output-json', required=True,
                        help="Path to output JSON summary")
    parser.add_argument('--output-html', required=True,
                        help="Path to output HTML report")
    
    # Parameters
    parser.add_argument('--segment', default="",
                        help="Segment name (for segmented genomes)")
    parser.add_argument('--pathogen', default="pathogen",
                        help="Pathogen name")
    parser.add_argument('--min-qc-score', type=float, default=80,
                        help="Minimum QC score threshold")
    
    return parser.parse_args()

def count_sequences(fasta_file):
    """Count number of sequences in FASTA file"""
    try:
        count = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        return count
    except Exception as e:
        logger.error(f"Error counting sequences in {fasta_file}: {e}")
        return 0

def extract_filter_reasons(log_file):
    """Extract filtering reasons from augur filter log"""
    filter_reasons = {}
    
    try:
        with open(log_file, "r") as f:
            log_text = f.read()
            
        # Example pattern for extracting info from augur filter log
        excluded_pattern = r"Excluding sequences: ([0-9]+) sequences excluded by filters"
        excluded_match = re.search(excluded_pattern, log_text)
        if excluded_match:
            filter_reasons["total_excluded"] = int(excluded_match.group(1))
        
        # Extract sequence length filter results
        length_pattern = r"filter: excluding ([0-9]+) sequences due to length \(min_length: ([0-9]+)\)"
        length_match = re.search(length_pattern, log_text)
        if length_match:
            filter_reasons["length"] = {
                "excluded": int(length_match.group(1)),
                "threshold": int(length_match.group(2))
            }
        
        # Extract ambiguous nucleotide filter results
        ambig_pattern = r"filter: excluding ([0-9]+) sequences due to excess ambiguous nucleotides \(max allowed: ([0-9]+)\)"
        ambig_match = re.search(ambig_pattern, log_text)
        if ambig_match:
            filter_reasons["ambiguous"] = {
                "excluded": int(ambig_match.group(1)),
                "threshold": int(ambig_match.group(2))
            }
            
    except Exception as e:
        logger.error(f"Error extracting filter reasons from {log_file}: {e}")
    
    return filter_reasons

def extract_nextclade_stats(json_file):
    """Extract statistics from Nextclade JSON output"""
    nextclade_stats = {
        "clade_distribution": {},
        "qc_scores": {},
        "mutation_counts": {},
        "missing_data": {}
    }
    
    try:
        with open(json_file, "r") as f:
            data = json.load(f)
        
        for result in data.get("results", []):
            # Clade information (if available)
            clade = result.get("clade", "Unknown")
            if clade not in nextclade_stats["clade_distribution"]:
                nextclade_stats["clade_distribution"][clade] = 0
            nextclade_stats["clade_distribution"][clade] += 1
            
            # QC score ranges
            qc_score = result.get("qc", {}).get("overallScore", 0)
            score_range = f"{int(qc_score // 10) * 10}-{int(qc_score // 10) * 10 + 9}"
            if score_range not in nextclade_stats["qc_scores"]:
                nextclade_stats["qc_scores"][score_range] = 0
            nextclade_stats["qc_scores"][score_range] += 1
            
            # Mutation counts
            substCounts = len(result.get("substitutions", []))
            if substCounts not in nextclade_stats["mutation_counts"]:
                nextclade_stats["mutation_counts"][substCounts] = 0
            nextclade_stats["mutation_counts"][substCounts] += 1
            
            # Missing data
            missing_data_pct = result.get("qc", {}).get("missingData", {}).get("score", 0)
            missing_range = f"{int(missing_data_pct * 10) / 10:.1f}-{int(missing_data_pct * 10) / 10 + 0.1:.1f}"
            if missing_range not in nextclade_stats["missing_data"]:
                nextclade_stats["missing_data"][missing_range] = 0
            nextclade_stats["missing_data"][missing_range] += 1
            
    except Exception as e:
        logger.error(f"Error extracting Nextclade stats from {json_file}: {e}")
    
    return nextclade_stats

def generate_html_report(qc_data, output_file):
    """Generate an HTML report from the QC data"""
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>RVF Sequence QC Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1, h2, h3 {{ color: #2c3e50; }}
        .summary {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        .section {{ margin-bottom: 30px; }}
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
    <h1>Sequence QC Report: {qc_data['pathogen_name']}</h1>
    <p>Generated on: {qc_data['report_date']}</p>
    <p>Segment: {qc_data['segment'] if qc_data['segment'] else 'All'}</p>
    
    <div class="summary">
        <h2>Summary</h2>
        <div class="metric">
            <span class="metric-name">Initial sequences:</span> {qc_data['sequence_counts']['raw']}
        </div>
        <div class="metric">
            <span class="metric-name">After basic filtering:</span> {qc_data['sequence_counts']['filtered']} 
            ({qc_data['sequence_counts']['filtered'] / qc_data['sequence_counts']['raw'] * 100:.1f}%)
        </div>
        <div class="metric">
            <span class="metric-name">After Nextclade QC:</span> {qc_data['sequence_counts']['nextclade_passed']} 
            ({qc_data['sequence_counts']['nextclade_passed'] / qc_data['sequence_counts']['raw'] * 100:.1f}%)
        </div>
    </div>
    
    <div class="section">
        <h2>Basic Filtering Results</h2>
        <table>
            <tr>
                <th>Filter Type</th>
                <th>Sequences Excluded</th>
                <th>Threshold</th>
            </tr>
    """
    
    # Add rows for each filter type
    for filter_type, details in qc_data.get('filter_reasons', {}).items():
        if isinstance(details, dict) and 'excluded' in details:
            html_content += f"""
            <tr>
                <td>{filter_type.capitalize()}</td>
                <td>{details['excluded']}</td>
                <td>{details.get('threshold', 'N/A')}</td>
            </tr>
            """
    
    html_content += f"""
        </table>
    </div>
    
    <div class="section">
        <h2>Nextclade QC Results</h2>
        
        <h3>QC Score Distribution</h3>
        <table>
            <tr>
                <th>Score Range</th>
                <th>Count</th>
            </tr>
    """
    
    # Add rows for QC score distribution
    for score_range, count in sorted(qc_data.get('nextclade_stats', {}).get('qc_scores', {}).items()):
        html_content += f"""
            <tr>
                <td>{score_range}</td>
                <td>{count}</td>
            </tr>
        """
    
    html_content += f"""
        </table>
        
        <h3>Clade Distribution</h3>
        <table>
            <tr>
                <th>Clade</th>
                <th>Count</th>
            </tr>
    """
    
    # Add rows for clade distribution
    for clade, count in qc_data.get('nextclade_stats', {}).get('clade_distribution', {}).items():
        html_content += f"""
            <tr>
                <td>{clade}</td>
                <td>{count}</td>
            </tr>
        """
    
    html_content += f"""
        </table>
    </div>
    
    <div class="section">
        <h2>Data Quality Metrics</h2>
        
        <h3>Missing Data Distribution</h3>
        <table>
            <tr>
                <th>Missing Data Range</th>
                <th>Count</th>
            </tr>
    """
    
    # Add rows for missing data distribution
    for missing_range, count in sorted(qc_data.get('nextclade_stats', {}).get('missing_data', {}).items()):
        html_content += f"""
            <tr>
                <td>{missing_range}</td>
                <td>{count}</td>
            </tr>
        """
    
    html_content += f"""
        </table>
        
        <h3>Mutation Count Distribution</h3>
        <table>
            <tr>
                <th>Number of Mutations</th>
                <th>Count</th>
            </tr>
    """
    
    # Add rows for mutation count distribution
    for mut_count, count in sorted(qc_data.get('nextclade_stats', {}).get('mutation_counts', {}).items()):
        html_content += f"""
            <tr>
                <td>{mut_count}</td>
                <td>{count}</td>
            </tr>
        """
    
    html_content += """
        </table>
    </div>
    
</body>
</html>
    """
    
    try:
        with open(output_file, 'w') as f:
            f.write(html_content)
        logger.info(f"HTML report written to {output_file}")
    except Exception as e:
        logger.error(f"Error writing HTML report: {e}")

def main():
    args = parse_args()
    
    try:
        # Count sequences at different stages
        raw_count = count_sequences(args.raw_sequences)
        filtered_count = count_sequences(args.filtered_sequences)
        nextclade_passed_count = count_sequences(args.nextclade_passed)
        
        # Extract filter reasons from logs
        filter_reasons = extract_filter_reasons(args.filter_log)
        
        # Extract Nextclade statistics
        nextclade_stats = extract_nextclade_stats(args.nextclade_json)
        
        # Calculate metadata statistics
        raw_metadata = pd.read_csv(args.raw_metadata, sep='\t')
        filtered_metadata = pd.read_csv(args.filtered_metadata, sep='\t')
        
        # Compile QC data for reporting
        qc_data = {
            "pathogen_name": args.pathogen,
            "segment": args.segment,
            "report_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "sequence_counts": {
                "raw": raw_count,
                "filtered": filtered_count,
                "nextclade_passed": nextclade_passed_count
            },
            "metadata_counts": {
                "raw": len(raw_metadata),
                "filtered": len(filtered_metadata)
            },
            "filter_reasons": filter_reasons,
            "nextclade_stats": nextclade_stats,
            "qc_thresholds": {
                "min_qc_score": args.min_qc_score
            }
        }
        
        # Write JSON summary
        with open(args.output_json, 'w') as f:
            json.dump(qc_data, f, indent=2)
        logger.info(f"QC summary JSON written to {args.output_json}")
        
        # Generate HTML report
        generate_html_report(qc_data, args.output_html)
        
    except Exception as e:
        logger.error(f"Error generating QC report: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()