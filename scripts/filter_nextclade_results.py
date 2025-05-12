#!/usr/bin/env python3
"""
Filters sequences based on Nextclade QC results.
"""

import os
import sys
import json
import argparse
import logging
from Bio import SeqIO

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter sequences based on Nextclade QC results")
    parser.add_argument('--input-json', required=True,
                        help="Path to Nextclade JSON output")
    parser.add_argument('--input-fasta', required=True,
                        help="Path to input sequences FASTA file")
    parser.add_argument('--output-fasta', required=True,
                        help="Path to write filtered sequences")
    parser.add_argument('--min-qc-score', type=float, default=80,
                        help="Minimum QC score to pass filtering")
    parser.add_argument('--segment', default="",
                        help="Segment name (for segmented genomes)")
    return parser.parse_args()

def filter_by_qc_score(nextclade_data, min_qc_score):
    """
    Extract sequence IDs that pass the QC threshold
    """
    passed_ids = set()
    failed_reasons = {}
    
    for result in nextclade_data.get("results", []):
        seq_id = result.get("seqName")
        qc_score = result.get("qc", {}).get("overallScore", 0)
        qc_status = result.get("qc", {}).get("overallStatus", "bad")
        
        if qc_score >= min_qc_score:
            passed_ids.add(seq_id)
        else:
            # Track failure reasons for reporting
            failed_reasons[seq_id] = {
                "score": qc_score,
                "status": qc_status,
                "issues": [
                    issue["tagValue"] 
                    for issue in result.get("qc", {}).get("privateMutations", {}).get("issues", []) +
                             result.get("qc", {}).get("missingData", {}).get("issues", []) +
                             result.get("qc", {}).get("snpClusters", {}).get("issues", [])
                ]
            }
    
    return passed_ids, failed_reasons

def main():
    args = parse_args()
    
    try:
        # Load Nextclade results
        with open(args.input_json) as f:
            nextclade_data = json.load(f)
            
        # Get sequences that pass QC
        passed_ids, failed_reasons = filter_by_qc_score(nextclade_data, args.min_qc_score)
        
        # Log QC statistics
        total_seqs = len(nextclade_data.get("results", []))
        passed_count = len(passed_ids)
        failed_count = total_seqs - passed_count
        
        logger.info(f"Nextclade QC summary for {args.segment if args.segment else 'all'} segment:")
        logger.info(f"  Total sequences: {total_seqs}")
        logger.info(f"  Passed QC (score >= {args.min_qc_score}): {passed_count}")
        logger.info(f"  Failed QC: {failed_count}")
        
        # Write passed sequences to output FASTA
        with open(args.output_fasta, "w") as out_f:
            for record in SeqIO.parse(args.input_fasta, "fasta"):
                if record.id in passed_ids:
                    SeqIO.write(record, out_f, "fasta")
        
        # Log details about failed sequences
        if failed_count > 0:
            logger.info("Failed sequences and reasons:")
            for seq_id, details in failed_reasons.items():
                issues = ", ".join(details["issues"]) if details["issues"] else "Low overall score"
                logger.info(f"  {seq_id}: Score={details['score']}, Issues: {issues}")
        
        logger.info(f"Filtered sequences written to {args.output_fasta}")
        
    except Exception as e:
        logger.error(f"Error in Nextclade filtering: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()