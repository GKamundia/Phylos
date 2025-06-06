#!/usr/bin/env python3
"""
Very simple direct filter for RVF sequences based on Nuc_Completeness
"""

import os
import sys
import pandas as pd
import traceback
from Bio import SeqIO

# Define file paths
metadata_path = "data/metadata/rvf_metadata_with_segments.tsv"
sequences_path = "data/sequences/raw/rvf_sequences.fasta"
output_metadata_path = "results/filtered/rvf_metadata.tsv"
output_sequences_path = "results/filtered/rvf_filtered.fasta"
log_path = "logs/filter_rvf.log"

try:
    # Create output directory
    os.makedirs("results/filtered", exist_ok=True)

    # Read metadata
    metadata = pd.read_csv(metadata_path, sep='\t', encoding='utf-8')
    total_count = len(metadata)
    
    # Filter metadata to keep only complete sequences
    filtered_metadata = metadata[metadata['Nuc_Completeness'] == 'complete']
    kept_count = len(filtered_metadata)
    dropped_count = total_count - kept_count
    
    # Get the list of sequence IDs to keep
    keep_ids = set(filtered_metadata['strain'])
    
    # Write filtered metadata
    filtered_metadata.to_csv(output_metadata_path, sep='\t', index=False, encoding='utf-8')
    
    # Filter sequences
    sequence_count = 0
    with open(sequences_path, "r", encoding='utf-8') as input_handle, open(output_sequences_path, "w", encoding='utf-8') as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if record.id in keep_ids:
                SeqIO.write(record, output_handle, "fasta")
                sequence_count += 1
    
    # Write log
    with open(log_path, 'w', encoding='utf-8') as log_file:
        log_file.write("Filtering for complete sequences across all segments\n")
        log_file.write(f"Total sequences before filtering: {total_count}\n")
        log_file.write(f"Sequences retained after filtering: {kept_count}\n")
        log_file.write(f"Sequences dropped due to 'Nuc_Completeness != complete': {dropped_count}\n")
        log_file.write(f"Filtering complete. Written {sequence_count} sequences to {output_sequences_path}\n")
        
        # Segment distribution
        if 'Segment' in filtered_metadata.columns:
            segment_counts = filtered_metadata['Segment'].value_counts().to_dict()
            log_file.write("Segment distribution in filtered data:\n")
            for segment, count in segment_counts.items():
                log_file.write(f"  {segment}: {count}\n")
    
    # Write a success marker file
    with open("filter_success.txt", "w", encoding='utf-8') as success_file:
        success_file.write("Filtering completed successfully\n")
        success_file.write(f"Total sequences before: {total_count}\n")
        success_file.write(f"Sequences kept: {kept_count}\n")
        success_file.write(f"Sequences filtered out: {dropped_count}\n")

except Exception as e:
    error_message = f"Error during filtering: {str(e)}\n{traceback.format_exc()}"
    print(error_message, file=sys.stderr)
    
    # Try to write to the log file
    try:
        with open(log_path, 'w', encoding='utf-8') as log_file:
            log_file.write(error_message)
    except:
        pass
    
    # Write to an error file
    with open("filter_error.txt", "w", encoding='utf-8') as error_file:
        error_file.write(error_message)
