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

# Dictionary to map output paths to segments
output_seq_map = {segment: path for segment, path in zip(segments, output_sequences)}
output_meta_map = {segment: path for segment, path in zip(segments, output_metadata)}

# Load metadata
logger.info(f"Loading metadata from {metadata_file}")
metadata = pd.read_csv(metadata_file, sep='\t')

# Split metadata by segment
segment_metadata = {}
for segment in segments:
    # Make segment matching case-insensitive
    segment_lower = segment.lower()
    segment_metadata[segment] = metadata[metadata['segment'] == segment_lower].copy()
    os.makedirs(os.path.dirname(output_meta_map[segment]), exist_ok=True)
    segment_metadata[segment].to_csv(output_meta_map[segment], sep='\t', index=False)
    logger.info(f"Wrote {len(segment_metadata[segment])} records to {output_meta_map[segment]}")

# Create segment-specific sequence files
for segment in segments:
    # Use accession IDs instead of strain names for matching sequence IDs
    segment_ids = set(segment_metadata[segment]['accession'].tolist())
    segment_seqs = []
    
    # Create output directory
    os.makedirs(os.path.dirname(output_seq_map[segment]), exist_ok=True)
    
    logger.info(f"Processing sequences for segment {segment}")
    
    # Process sequences
    seq_count = 0
    with open(output_seq_map[segment], 'w') as outfile:
        for record in SeqIO.parse(sequences_file, "fasta"):
            strain = record.id
            if strain in segment_ids:
                SeqIO.write(record, outfile, "fasta")
                seq_count += 1
    
    logger.info(f"Wrote {seq_count} sequences to {output_seq_map[segment]}")

logger.info("Segment splitting completed successfully")