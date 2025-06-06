#!/usr/bin/env python3
"""
Debug script to verify we can read and write files
"""

import sys
import os
import argparse
import pandas as pd
from Bio import SeqIO

def main():
    print("DEBUG: Starting debug filter script")
    
    try:
        # Print Python version and working directory
        print(f"Python version: {sys.version}")
        print(f"Current working directory: {os.getcwd()}")
        
        # Test file reading
        metadata_path = "data/metadata/rvf_metadata_with_segments.tsv"
        sequences_path = "data/sequences/raw/rvf_sequences.fasta"
        
        # Check if files exist
        print(f"Metadata file exists: {os.path.exists(metadata_path)}")
        print(f"Sequences file exists: {os.path.exists(sequences_path)}")
        
        # Create output directory
        output_dir = "results/filtered"
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
        
        # Test metadata reading
        if os.path.exists(metadata_path):
            print(f"Reading metadata from {metadata_path}")
            metadata = pd.read_csv(metadata_path, sep='\t')
            print(f"Metadata shape: {metadata.shape}")
            print(f"Metadata columns: {list(metadata.columns)}")
            
            # Check for Nuc_Completeness column
            if 'Nuc_Completeness' in metadata.columns:
                print("Nuc_Completeness column found in metadata")
                complete_count = len(metadata[metadata['Nuc_Completeness'] == 'complete'])
                print(f"Complete sequences: {complete_count}")
                print(f"Partial sequences: {len(metadata) - complete_count}")
                
                # Write filtered metadata
                filtered_metadata = metadata[metadata['Nuc_Completeness'] == 'complete']
                output_metadata_path = os.path.join(output_dir, "rvf_metadata.tsv")
                print(f"Writing filtered metadata to {output_metadata_path}")
                filtered_metadata.to_csv(output_metadata_path, sep='\t', index=False)
                
                # Test sequence reading
                if os.path.exists(sequences_path):
                    print(f"Reading sequences from {sequences_path}")
                    sequence_count = 0
                    output_sequences_path = os.path.join(output_dir, "rvf_filtered.fasta")
                    
                    # Get the list of sequence IDs to keep
                    keep_ids = set(filtered_metadata['strain'])
                    print(f"Keeping {len(keep_ids)} unique strain IDs")
                    
                    with open(sequences_path, "r") as input_handle, open(output_sequences_path, "w") as output_handle:
                        for record in SeqIO.parse(input_handle, "fasta"):
                            if record.id in keep_ids:
                                SeqIO.write(record, output_handle, "fasta")
                                sequence_count += 1
                                if sequence_count % 100 == 0:
                                    print(f"Processed {sequence_count} sequences")
                    
                    print(f"Filtering complete. Written {sequence_count} sequences to {output_sequences_path}")
                else:
                    print(f"ERROR: Could not read sequences file at {sequences_path}")
            else:
                print(f"ERROR: Nuc_Completeness column not found in metadata!")
                print(f"Available columns: {list(metadata.columns)}")
        else:
            print(f"ERROR: Could not read metadata file at {metadata_path}")
        
        # Write a simple log file
        with open("logs/debug_filter.log", "w") as log_file:
            log_file.write("Debug filter script completed successfully\n")
        
    except Exception as e:
        print(f"ERROR during debugging: {str(e)}")
        with open("logs/debug_filter.log", "w") as log_file:
            log_file.write(f"Error during debugging: {str(e)}\n")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
