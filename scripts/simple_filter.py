#!/usr/bin/env python3
"""
Simple direct filter for RVF sequences based on Nuc_Completeness
"""

import os
import pandas as pd
from Bio import SeqIO

# Define file paths
metadata_path = "data/metadata/rvf_metadata_with_segments.tsv"
sequences_path = "data/sequences/raw/rvf_sequences.fasta"
output_metadata_path = "results/filtered/rvf_metadata.tsv"
output_sequences_path = "results/filtered/rvf_filtered.fasta"
log_path = "logs/filter_rvf.log"

# Create output directory
os.makedirs("results/filtered", exist_ok=True)

# Open log file
with open(log_path, 'w') as log_file:
    log_file.write("Filtering for complete sequences across all segments\n")
    
    try:
        # Read metadata
        metadata = pd.read_csv(metadata_path, sep='\t')
        total_count = len(metadata)
        log_file.write(f"Total sequences before filtering: {total_count}\n")
        
        # Filter metadata to keep only complete sequences
        filtered_metadata = metadata[metadata['Nuc_Completeness'] == 'complete']
        kept_count = len(filtered_metadata)
        dropped_count = total_count - kept_count
        
        log_file.write(f"Sequences retained after filtering: {kept_count}\n")
        log_file.write(f"Sequences dropped due to 'Nuc_Completeness != complete': {dropped_count}\n")
        
        # Get the list of sequence IDs to keep
        keep_ids = set(filtered_metadata['strain'])
        
        # Write filtered metadata
        filtered_metadata.to_csv(output_metadata_path, sep='\t', index=False)
        
        # Filter sequences
        sequence_count = 0
        with open(sequences_path, "r") as input_handle, open(output_sequences_path, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                if record.id in keep_ids:
                    SeqIO.write(record, output_handle, "fasta")
                    sequence_count += 1
        
        log_file.write(f"Filtering complete. Written {sequence_count} sequences to {output_sequences_path}\n")
        
        # Segment distribution
        if 'Segment' in filtered_metadata.columns:
            segment_counts = filtered_metadata['Segment'].value_counts().to_dict()
            log_file.write("Segment distribution in filtered data:\n")
            for segment, count in segment_counts.items():
                log_file.write(f"  {segment}: {count}\n")
        
    except Exception as e:
        log_file.write(f"Error during filtering: {str(e)}\n")
