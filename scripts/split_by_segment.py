#!/usr/bin/env python3
"""
Split sequences by segment for multi-segment pathogen analysis
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
    
    # Process sequences
    seq_count = 0
    with open(output_seq_map[segment], 'w') as outfile:
        for record in sequences:
            accession = record.id
            if accession in segment_accessions:
                SeqIO.write(record, outfile, "fasta")
                seq_count += 1
    
    logger.info(f"Wrote {seq_count} sequences to {output_seq_map[segment]}")
    
    # Warn if there's a mismatch
    if seq_count != len(segment_accessions):
        missing_count = len(segment_accessions) - seq_count
        logger.warning(f"Segment {segment}: Found {seq_count} sequences but expected {len(segment_accessions)} (missing: {missing_count})")

# Create done file if it's specified in output
try:
    done_file = snakemake.output.done
    with open(done_file, 'w') as f:
        f.write("Segment splitting completed successfully\n")
        f.write(f"Processed segments: {', '.join(segments)}\n")
        for segment in segments:
            seg_count = len(segment_metadata[segment])
            f.write(f"{segment}: {seg_count} sequences\n")
    logger.info(f"Created done file: {done_file}")
except AttributeError:
    # No done file specified
    logger.info("No done file specified in output")
    pass

logger.info("Segment splitting completed successfully")