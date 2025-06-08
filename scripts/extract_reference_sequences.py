#!/usr/bin/env python3
"""
Extract NCBI RefSeq reference sequences for RVF segments from raw data.
These will be used as alignment references for each segment.
"""

import argparse
from Bio import SeqIO
import os

def extract_reference_sequences(raw_sequences_file, output_dir):
    """
    Extract NCBI RefSeq reference sequences for RVF segments.
    
    Reference sequences:
    - NC_014395.1: S segment (1690 bp)
    - NC_014396.1: M segment (3885 bp) 
    - NC_014397.1: L segment (6404 bp)
    """
    
    # Define the RefSeq accessions for each segment
    refseq_ids = {
        'NC_014395.1': 'S',  # S segment
        'NC_014396.1': 'M',  # M segment  
        'NC_014397.1': 'L'   # L segment
    }
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Track found sequences
    found_sequences = {}
    
    # Parse the raw sequences file
    print(f"Extracting reference sequences from {raw_sequences_file}")
    
    with open(raw_sequences_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in refseq_ids:
                segment = refseq_ids[record.id]
                found_sequences[segment] = record
                print(f"Found {segment} segment reference: {record.id} ({len(record.seq)} bp)")
    
    # Write reference files for each segment
    for segment in ['L', 'M', 'S']:
        output_file = os.path.join(output_dir, f"reference_{segment.lower()}.fasta")
        
        if segment in found_sequences:
            with open(output_file, 'w') as output_handle:
                SeqIO.write(found_sequences[segment], output_handle, 'fasta')
            print(f"Reference for segment {segment} written to {output_file}")
        else:
            print(f"WARNING: Reference for segment {segment} not found!")
    
    # Also create a combined reference file
    combined_output = os.path.join(output_dir, "reference_all_segments.fasta")
    with open(combined_output, 'w') as output_handle:
        for segment in ['L', 'M', 'S']:
            if segment in found_sequences:
                SeqIO.write(found_sequences[segment], output_handle, 'fasta')
    
    print(f"Combined reference file written to {combined_output}")
    
    return found_sequences

def main():
    parser = argparse.ArgumentParser(description='Extract NCBI RefSeq reference sequences for RVF')
    parser.add_argument('--sequences', required=True, help='Path to raw sequences FASTA file')
    parser.add_argument('--output-dir', required=True, help='Output directory for reference files')
    
    args = parser.parse_args()
    
    extract_reference_sequences(args.sequences, args.output_dir)

if __name__ == "__main__":
    main()
