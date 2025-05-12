#!/usr/bin/env python3
"""
Validates FASTA files for correctness and provides statistics
"""

import os
import sys
import argparse
import logging
from Bio import SeqIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Validate FASTA file format and content")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output validation report")
    return parser.parse_args()

def validate_fasta(fasta_file):
    """Validate a FASTA file and return statistics and any issues found"""
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        
        if not records:
            return False, {"error": "No sequences found in the FASTA file"}, {}
        
        # Calculate statistics
        stats = {
            "num_sequences": len(records),
            "min_length": min(len(record.seq) for record in records),
            "max_length": max(len(record.seq) for record in records),
            "avg_length": sum(len(record.seq) for record in records) / len(records),
            "unique_ids": len(set(record.id for record in records))
        }
        
        # Check for issues
        issues = []
        sequence_ids = set()
        
        for i, record in enumerate(records):
            # Check for duplicate IDs
            if record.id in sequence_ids:
                issues.append(f"Duplicate sequence ID found: {record.id}")
            sequence_ids.add(record.id)
            
            # Check for empty sequences
            if len(record.seq) == 0:
                issues.append(f"Empty sequence found at index {i}, ID: {record.id}")
            
            # Check for non-nucleotide characters
            valid_chars = set('ATGCNatgcn-')
            invalid_chars = set(str(record.seq)) - valid_chars
            if invalid_chars:
                issues.append(f"Invalid nucleotide characters found in sequence {record.id}: {', '.join(invalid_chars)}")
        
        # Check for duplicate IDs
        if stats["unique_ids"] != stats["num_sequences"]:
            issues.append(f"Found {stats['num_sequences'] - stats['unique_ids']} duplicate sequence IDs")
        
        return len(issues) == 0, stats, issues
        
    except Exception as e:
        return False, {"error": str(e)}, {}

def main():
    args = parse_args()
    
    logger.info(f"Validating FASTA file: {args.input}")
    
    valid, stats, issues = validate_fasta(args.input)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Write validation report
    with open(args.output, "w") as f:
        f.write("FASTA Validation Report\n")
        f.write("=====================\n\n")
        
        f.write(f"File: {args.input}\n")
        f.write(f"Valid: {'Yes' if valid else 'No'}\n\n")
        
        f.write("Statistics:\n")
        for key, value in stats.items():
            if key != "error":  # Skip error if it exists
                f.write(f"  {key}: {value}\n")
        
        if issues:
            f.write("\nIssues:\n")
            for issue in issues:
                f.write(f"  * {issue}\n")
        else:
            f.write("\nNo issues found.\n")
    
    if valid:
        logger.info(f"FASTA validation successful. Found {stats['num_sequences']} sequences.")
        return 0
    else:
        if "error" in stats:
            logger.error(f"FASTA validation failed: {stats['error']}")
        else:
            logger.error(f"FASTA validation failed with {len(issues)} issues.")
        return 1

if __name__ == "__main__":
    sys.exit(main())