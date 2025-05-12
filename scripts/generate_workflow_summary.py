#!/usr/bin/env python3
"""
Generates a summary of the Nextstrain workflow run
"""

import os
import sys
import json
import time
from datetime import datetime
from pathlib import Path
from Bio import SeqIO

# Access Snakemake variables
pathogen = snakemake.params.pathogen
pathogen_name = snakemake.params.pathogen_name
segment_mode = snakemake.params.segment_mode
segments = snakemake.params.segments
summary_output = snakemake.output.summary

def count_sequences(fasta_file):
    """Count sequences in a FASTA file if it exists"""
    try:
        if os.path.exists(fasta_file):
            return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        else:
            return 0
    except Exception:
        return 0

def get_file_info(filepath):
    """Get file size and modification time for a file"""
    try:
        if os.path.exists(filepath):
            stats = os.stat(filepath)
            return {
                "size_bytes": stats.st_size,
                "size_human": format_file_size(stats.st_size),
                "modified": datetime.fromtimestamp(stats.st_mtime).isoformat()
            }
        else:
            return {
                "size_bytes": 0,
                "size_human": "0 B",
                "modified": None
            }
    except Exception:
        return {
            "size_bytes": 0, 
            "size_human": "0 B", 
            "modified": None
        }

def format_file_size(size_bytes):
    """Format file size in human-readable format"""
    for unit in ["B", "KB", "MB", "GB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.2f} TB"

def main():
    # Initialize summary dictionary
    summary = {
        "pathogen": pathogen,
        "pathogen_name": pathogen_name,
        "segment_mode": segment_mode,
        "segments": segments,
        "timestamp": datetime.now().isoformat(),
        "sequence_counts": {},
        "files": {}
    }
    
    # Collect sequence counts at different stages
    if segment_mode == "single":
        raw_fasta = f"data/sequences/raw/{pathogen}_sequences.fasta"
        filtered_fasta = f"results/filtered/{pathogen}_filtered.fasta"
        aligned_fasta = f"results/aligned/{pathogen}_aligned.fasta"
        nextclade_passed = f"results/nextclade/{pathogen}_passed.fasta"
        
        summary["sequence_counts"] = {
            "raw": count_sequences(raw_fasta),
            "filtered": count_sequences(filtered_fasta),
            "aligned": count_sequences(aligned_fasta),
            "qc_passed": count_sequences(nextclade_passed),
        }
        
        # Collect file information
        output_json = f"results/auspice/{pathogen}.json"
        qc_report = f"results/qc_reports/{pathogen}_qc_report.html"
        
        summary["files"] = {
            "input_fasta": get_file_info(raw_fasta),
            "output_json": get_file_info(output_json),
            "qc_report": get_file_info(qc_report)
        }
    else:
        # Multi-segment mode
        raw_fasta = f"data/sequences/raw/{pathogen}_sequences.fasta"
        summary["sequence_counts"]["raw"] = count_sequences(raw_fasta)
        
        segment_counts = {}
        for segment in segments:
            segment_filtered = f"results/segments/{segment}/filtered/{pathogen}_{segment}_filtered.fasta"
            segment_aligned = f"results/segments/{segment}/aligned/{pathogen}_{segment}_aligned.fasta"
            segment_qc = f"results/segments/{segment}/nextclade/{pathogen}_{segment}_passed.fasta"
            
            segment_counts[segment] = {
                "filtered": count_sequences(segment_filtered),
                "aligned": count_sequences(segment_aligned),
                "qc_passed": count_sequences(segment_qc)
            }
        
        summary["segment_counts"] = segment_counts
        
        # Collect file information for combined output
        combined_json = f"results/auspice/{pathogen}_combined.json"
        combined_qc = f"results/qc_reports/{pathogen}_combined_qc_report.html"
        
        summary["files"] = {
            "input_fasta": get_file_info(raw_fasta),
            "combined_json": get_file_info(combined_json),
            "combined_qc_report": get_file_info(combined_qc)
        }
        
        # Add segment-specific files
        for segment in segments:
            segment_json = f"results/segments/{segment}/auspice/{pathogen}_{segment}.json"
            segment_qc = f"results/segments/{segment}/qc_reports/{pathogen}_{segment}_qc_report.html"
            
            summary["files"][f"segment_{segment}_json"] = get_file_info(segment_json)
            summary["files"][f"segment_{segment}_qc"] = get_file_info(segment_qc)
    
    # Calculate total runtime (approximate based on file timestamps)
    raw_meta_file = f"data/metadata/raw/{pathogen}_metadata.tsv"
    final_output = summary["files"].get("output_json", summary["files"].get("combined_json", {}))
    
    try:
        if os.path.exists(raw_meta_file) and final_output.get("modified"):
            start_time = os.stat(raw_meta_file).st_mtime
            end_time = datetime.fromisoformat(final_output["modified"]).timestamp()
            runtime_seconds = end_time - start_time
            
            # Format runtime
            hours, remainder = divmod(runtime_seconds, 3600)
            minutes, seconds = divmod(remainder, 60)
            
            summary["runtime"] = {
                "seconds": runtime_seconds,
                "formatted": f"{int(hours)}h {int(minutes)}m {int(seconds)}s"
            }
    except Exception:
        summary["runtime"] = {
            "seconds": 0,
            "formatted": "unknown"
        }
    
    # Write summary to file
    os.makedirs(os.path.dirname(summary_output), exist_ok=True)
    with open(summary_output, "w") as f:
        json.dump(summary, f, indent=2)

if __name__ == "__main__":
    main()