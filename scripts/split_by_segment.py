#!/usr/bin/env python3
"""
Split sequences by segment for multi-segment pathogen analysis.
Only includes sequences that match exact reference le            for length, count in sorted_lengths[:5]:  # Show top 5 lengths
                status = "KEPT" if length == reference_length else "FILTERED"
                logger.info(f"    {length}bp: {count} sequences {status}")hs for better phylogenetic analysis.
"""

import os
import pandas as pd
from Bio import SeqIO
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Define exact reference lengths for filtering
REFERENCE_LENGTHS = {
    'L': 6404,  # L segment exact reference length  
    'M': 3885,  # M segment exact reference length
    'S': 1690   # S segment exact reference length
}

# Get input parameters from Snakemake (automatically available in script mode)
sequences_file = snakemake.input.sequences
metadata_file = snakemake.input.metadata
segments = snakemake.params.segments
output_sequences = snakemake.output.sequences
output_metadata = snakemake.output.metadata
log_file = snakemake.log[0]

# Configure logging to file
file_handler = logging.FileHandler(log_file)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(file_handler)

logger.info(f"Starting segment splitting process")
logger.info(f"Input sequences: {sequences_file}")
logger.info(f"Input metadata: {metadata_file}")
logger.info(f"Target segments: {segments}")

# Dictionary to map output paths to segments
output_seq_map = {segment: path for segment, path in zip(segments, output_sequences)}
output_meta_map = {segment: path for segment, path in zip(segments, output_metadata)}

# Load metadata
logger.info(f"Loading metadata from {metadata_file}")
metadata = pd.read_csv(metadata_file, sep='\t')
logger.info(f"Loaded {len(metadata)} metadata records")

# Check segment column
if 'segment' not in metadata.columns:
    raise ValueError("Metadata must contain a 'segment' column")

logger.info(f"Segment distribution: {metadata['segment'].value_counts().to_dict()}")

def filter_sequences_by_length(sequences, segment, reference_length):
    """
    Filter sequences to only include those matching the exact reference length.
    
    Args:
        sequences: List of Bio.SeqRecord objects
        segment: Segment name (L, M, or S)
        reference_length: Expected reference length for this segment
        
    Returns:
        tuple: (filtered_sequences, length_stats)
    """
    filtered_sequences = []
    length_counts = {}
    
    for record in sequences:
        seq_length = len(record.seq)
        length_counts[seq_length] = length_counts.get(seq_length, 0) + 1
        
        if seq_length == reference_length:
            filtered_sequences.append(record)
    
    length_stats = {
        'total_input': len(sequences),
        'exact_match': len(filtered_sequences),
        'filter_rate': len(filtered_sequences) / len(sequences) * 100 if sequences else 0,
        'length_distribution': length_counts
    }
    
    logger.info(f"Segment {segment} length filtering:")
    logger.info(f"  Input sequences: {length_stats['total_input']}")
    logger.info(f"  Exact reference length ({reference_length}bp): {length_stats['exact_match']}")
    logger.info(f"  Retention rate: {length_stats['filter_rate']:.1f}%")
    
    return filtered_sequences, length_stats

# Split metadata by segment
segment_metadata = {}
for segment in segments:
    # Make segment matching case-insensitive
    segment_lower = segment.lower()
    segment_data = metadata[metadata['segment'] == segment_lower].copy()
    segment_metadata[segment] = segment_data
    
    # Create output directory and save metadata
    os.makedirs(os.path.dirname(output_meta_map[segment]), exist_ok=True)
    segment_data.to_csv(output_meta_map[segment], sep='\t', index=False)
    logger.info(f"Wrote {len(segment_data)} metadata records for segment {segment} to {output_meta_map[segment]}")

# Load sequences
logger.info(f"Loading sequences from {sequences_file}")
sequences = list(SeqIO.parse(sequences_file, "fasta"))
logger.info(f"Loaded {len(sequences)} sequences")

# Create segment-specific sequence files
for segment in segments:
    # Use accession IDs for matching sequence IDs
    segment_accessions = set(segment_metadata[segment]['accession'].tolist())
    
    # Create output directory
    os.makedirs(os.path.dirname(output_seq_map[segment]), exist_ok=True)
    
    logger.info(f"Processing sequences for segment {segment}")
    logger.info(f"Looking for {len(segment_accessions)} accessions in segment {segment}")
    
    # Get segment-specific sequences first
    segment_sequences = []
    for record in sequences:
        accession = record.id
        if accession in segment_accessions:
            segment_sequences.append(record)
    
    logger.info(f"Found {len(segment_sequences)} sequences for segment {segment}")
    
    # Apply exact length filtering
    reference_length = REFERENCE_LENGTHS.get(segment.upper())
    if reference_length:
        filtered_sequences, length_stats = filter_sequences_by_length(
            segment_sequences, segment, reference_length
        )
        
        # Log filtering details
        if length_stats['total_input'] > length_stats['exact_match']:
            rejected_count = length_stats['total_input'] - length_stats['exact_match']
            logger.info(f"  Filtered out {rejected_count} sequences with non-reference lengths")
              # Log the most common non-reference lengths
            sorted_lengths = sorted(length_stats['length_distribution'].items(), 
                                  key=lambda x: x[1], reverse=True)
            logger.info(f"  Length distribution:")
            for length, count in sorted_lengths[:5]:  # Show top 5 lengths
                status = "KEPT" if length == reference_length else "FILTERED"
                logger.info(f"    {length}bp: {count} sequences {status}")
    else:
        logger.warning(f"No reference length defined for segment {segment}")
        filtered_sequences = segment_sequences
    
    # Write filtered sequences to file
    seq_count = 0
    with open(output_seq_map[segment], 'w') as outfile:
        for record in filtered_sequences:
            SeqIO.write(record, outfile, "fasta")
            seq_count += 1
    
    logger.info(f"Wrote {seq_count} length-filtered sequences to {output_seq_map[segment]}")
    
    # Update metadata to only include sequences that passed length filtering
    filtered_accessions = {record.id for record in filtered_sequences}
    original_meta_count = len(segment_metadata[segment])
    segment_metadata[segment] = segment_metadata[segment][
        segment_metadata[segment]['accession'].isin(filtered_accessions)
    ].copy()
    
    # Re-save the filtered metadata
    segment_metadata[segment].to_csv(output_meta_map[segment], sep='\t', index=False)
    filtered_meta_count = len(segment_metadata[segment])
    
    logger.info(f"Updated metadata: {original_meta_count} -> {filtered_meta_count} records")
    
    # Warn if there's a significant reduction
    if seq_count != len(segment_accessions):
        reduction = len(segment_accessions) - seq_count
        reduction_pct = (reduction / len(segment_accessions)) * 100 if segment_accessions else 0
        logger.warning(f"Segment {segment}: Reduced from {len(segment_accessions)} to {seq_count} sequences ({reduction_pct:.1f}% reduction due to length filtering)")

# Create done file if it's specified in output
try:
    done_file = snakemake.output.done
    with open(done_file, 'w') as f:
        f.write("Segment splitting with exact length filtering completed successfully\n")
        f.write(f"Processed segments: {', '.join(segments)}\n")
        f.write(f"Reference lengths: L={REFERENCE_LENGTHS['L']}bp, M={REFERENCE_LENGTHS['M']}bp, S={REFERENCE_LENGTHS['S']}bp\n")
        f.write("Final sequence counts (after length filtering):\n")
        for segment in segments:
            seg_count = len(segment_metadata[segment])
            f.write(f"  {segment}: {seg_count} sequences\n")
    logger.info(f"Created done file: {done_file}")
except AttributeError:
    # No done file specified
    logger.info("No done file specified in output")
    pass

logger.info("Segment splitting with exact length filtering completed successfully")