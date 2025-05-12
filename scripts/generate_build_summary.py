#!/usr/bin/env python3
"""
Generate a summary report for a completed Nextstrain build
"""

import os
import json
import yaml
from datetime import datetime

def generate_build_summary():
    """Generate a summary of the latest build run."""
    try:
        # Read master config to get active pathogen
        with open("config/master_config.yaml", "r") as f:
            master_config = yaml.safe_load(f)
        
        pathogen = master_config["active_pathogen"]
        pathogen_name = master_config["pathogens"][pathogen]["name"]
        
        # Define output paths
        summary_path = f"results/build_summary_{pathogen}_{datetime.now().strftime('%Y%m%d')}.md"
        qc_summary_path = f"results/qc_reports/{pathogen}_qc_summary.json"
        
        # Check if QC summary exists
        if not os.path.exists(qc_summary_path):
            print(f"Warning: QC summary file not found at {qc_summary_path}")
            qc_data = {"sequence_counts": {}}
        else:
            with open(qc_summary_path, "r") as f:
                qc_data = json.load(f)
        
        # Get sequence counts
        seq_counts = qc_data.get("sequence_counts", {})
        raw_count = seq_counts.get("raw", "N/A")
        filtered_count = seq_counts.get("filtered", "N/A")
        qc_passed = seq_counts.get("qc_passed", "N/A")
        
        # Generate summary markdown
        summary = [
            f"# Nextstrain {pathogen_name} Build Summary",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Build Information",
            f"- Pathogen: {pathogen_name} ({pathogen})",
            f"- Build Date: {datetime.now().strftime('%Y-%m-%d')}",
            "",
            "## Sequence Counts",
            f"- Raw sequences: {raw_count}",
            f"- Filtered sequences: {filtered_count}",
            f"- QC passed sequences: {qc_passed}",
            "",
            "## Results",
            f"- Auspice JSON: results/auspice/{pathogen}.json",
            f"- QC Report: results/qc_reports/{pathogen}_qc_report.html",
            "",
            "## Next Steps",
            "- Review the results in Auspice: `nextstrain view results/auspice/`",
            "- Check the QC report for any data quality issues"
        ]
        
        # Write summary to file
        with open(summary_path, "w") as f:
            f.write("\n".join(summary))
        
        print(f"Build summary generated at {summary_path}")
        
    except Exception as e:
        print(f"Error generating build summary: {e}")

if __name__ == "__main__":
    generate_build_summary()