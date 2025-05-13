#!/usr/bin/env python3
"""
Update segment information in metadata based on sequence length
"""

import pandas as pd
import argparse
import sys
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Update segment information in metadata based on sequence length")
    parser.add_argument("--input-metadata", required=True, help="Input metadata TSV file")
    parser.add_argument("--input-sequences", required=True, help="Input sequences FASTA file")
    parser.add_argument("--output-metadata", required=True, help="Output metadata TSV file with updated segment information")
    return parser.parse_args()

def infer_segment_from_length(length):
    """
    Infer segment based on sequence length:
    
    L segment: ~6.4kb (typically >6000bp)
    M segment: ~3.9kb (typically 3500-6000bp)
    S segment: ~1.7kb (typically <3500bp)
    """
    if length >= 6000:
        return "L"
    elif length >= 3500 and length < 6000:
        return "M"
    elif length >= 1000:  # Ensure it's a reasonably sized sequence
        return "S"
    else:
        return "unknown"  # Keep as unknown for very short sequences

def main():
    args = parse_args()
    
    print(f"Reading metadata from {args.input_metadata}")
    metadata = pd.read_csv(args.input_metadata, sep='\t')
    
    print(f"Reading sequences from {args.input_sequences}")
    seq_lengths = {}
    for record in SeqIO.parse(args.input_sequences, "fasta"):
        seq_lengths[record.id] = len(record.seq)
    
    print(f"Found {len(seq_lengths)} sequences")
    
    # Update segment information based on sequence length
    updated_count = 0
    for idx, row in metadata.iterrows():
        strain = row['strain']
        
        # If we have the strain in our sequences, use that length
        if strain in seq_lengths:
            length = seq_lengths[strain]
            segment = infer_segment_from_length(length)
            metadata.at[idx, 'segment'] = segment
            updated_count += 1
        # Otherwise, try to use the length from metadata
        elif not pd.isna(row['length']) and row['length'] > 0:
            length = row['length']
            segment = infer_segment_from_length(length)
            metadata.at[idx, 'segment'] = segment
            updated_count += 1
    
    # Write updated metadata
    metadata.to_csv(args.output_metadata, sep='\t', index=False)
    
    # Print summary
    l_count = sum(metadata['segment'] == 'L')
    m_count = sum(metadata['segment'] == 'M')
    s_count = sum(metadata['segment'] == 'S')
    unknown_count = sum(metadata['segment'] == 'unknown')
    
    print(f"Updated segment information for {updated_count} sequences")
    print(f"L segment: {l_count}")
    print(f"M segment: {m_count}")
    print(f"S segment: {s_count}")
    print(f"Unknown: {unknown_count}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())