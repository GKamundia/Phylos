#!/usr/bin/env python3
"""
Standalone wrapper for split_by_segment.py to run outside Snakemake.
"""
import argparse
import os
import sys
from Bio import SeqIO
import pandas as pd
import logging

# --- Helper functions and constants (inlined for standalone use) ---
REFERENCE_LENGTHS = {
    'L': 6404,  # L segment exact reference length  
    'M': 3885,  # M segment exact reference length
    'S': 1690   # S segment exact reference length
}

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

# Set up argument Parser
parser = argparse.ArgumentParser(description="Split sequences and metadata by segment with exact length filtering.")
parser.add_argument('--sequences', required=True, help='Input FASTA file')
parser.add_argument('--metadata', required=True, help='Input metadata TSV file')
parser.add_argument('--segments', nargs='+', required=True, help='Segments to process (e.g. L M S)')
parser.add_argument('--output-sequences', nargs='+', required=True, help='Output FASTA files for each segment (order must match --segments)')
parser.add_argument('--output-metadata', nargs='+', required=True, help='Output metadata TSV files for each segment (order must match --segments)')
parser.add_argument('--log', required=True, help='Log file path')
args = parser.parse_args()

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
file_handler = logging.FileHandler(args.log)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(file_handler)

segments = args.segments
output_seq_map = {segment: path for segment, path in zip(segments, args.output_sequences)}
output_meta_map = {segment: path for segment, path in zip(segments, args.output_metadata)}

logger.info(f"Starting segment splitting process (standalone)")
logger.info(f"Input sequences: {args.sequences}")
logger.info(f"Input metadata: {args.metadata}")
logger.info(f"Target segments: {segments}")

# Load metadata
logger.info(f"Loading metadata from {args.metadata}")
metadata = pd.read_csv(args.metadata, sep='\t')
logger.info(f"Loaded {len(metadata)} metadata records")

if 'segment' not in metadata.columns:
    logger.error("Metadata must contain a 'segment' column")
    sys.exit(1)

logger.info(f"Segment distribution: {metadata['segment'].value_counts().to_dict()}")

# Split metadata by segment
segment_metadata = {}
for segment in segments:
    segment_lower = segment.lower()
    segment_data = metadata[metadata['segment'] == segment_lower].copy()
    segment_metadata[segment] = segment_data
    os.makedirs(os.path.dirname(output_meta_map[segment]), exist_ok=True)
    segment_data.to_csv(output_meta_map[segment], sep='\t', index=False)
    logger.info(f"Wrote {len(segment_data)} metadata records for segment {segment} to {output_meta_map[segment]}")

# Load sequences
logger.info(f"Loading sequences from {args.sequences}")
sequences = list(SeqIO.parse(args.sequences, "fasta"))
logger.info(f"Loaded {len(sequences)} sequences")

for segment in segments:
    segment_accessions = set(segment_metadata[segment]['accession'].tolist())
    os.makedirs(os.path.dirname(output_seq_map[segment]), exist_ok=True)
    logger.info(f"Processing sequences for segment {segment}")
    logger.info(f"Looking for {len(segment_accessions)} accessions in segment {segment}")
    segment_sequences = [record for record in sequences if record.id in segment_accessions]
    logger.info(f"Found {len(segment_sequences)} sequences for segment {segment}")
    reference_length = REFERENCE_LENGTHS.get(segment.upper())
    if reference_length:
        filtered_sequences, length_stats = filter_sequences_by_length(segment_sequences, segment, reference_length)
        if length_stats['total_input'] > length_stats['exact_match']:
            rejected_count = length_stats['total_input'] - length_stats['exact_match']
            logger.info(f"  Filtered out {rejected_count} sequences with non-reference lengths")
            sorted_lengths = sorted(length_stats['length_distribution'].items(), key=lambda x: x[1], reverse=True)
            logger.info(f"  Length distribution:")
            for length, count in sorted_lengths[:5]:
                status = "KEPT" if length == reference_length else "FILTERED"
                logger.info(f"    {length}bp: {count} sequences {status}")
    else:
        logger.warning(f"No reference length defined for segment {segment}")
        filtered_sequences = segment_sequences
    seq_count = 0
    with open(output_seq_map[segment], 'w') as outfile:
        for record in filtered_sequences:
            SeqIO.write(record, outfile, "fasta")
            seq_count += 1
    logger.info(f"Wrote {seq_count} length-filtered sequences to {output_seq_map[segment]}")
    filtered_accessions = {record.id for record in filtered_sequences}
    original_meta_count = len(segment_metadata[segment])
    segment_metadata[segment] = segment_metadata[segment][segment_metadata[segment]['accession'].isin(filtered_accessions)].copy()
    segment_metadata[segment].to_csv(output_meta_map[segment], sep='\t', index=False)
    filtered_meta_count = len(segment_metadata[segment])
    logger.info(f"Updated metadata: {original_meta_count} -> {filtered_meta_count} records")
    if seq_count != len(segment_accessions):
        reduction = len(segment_accessions) - seq_count
        reduction_pct = (reduction / len(segment_accessions)) * 100 if segment_accessions else 0
        logger.warning(f"Segment {segment}: Reduced from {len(segment_accessions)} to {seq_count} sequences ({reduction_pct:.1f}% reduction due to length filtering)")

logger.info("Segment splitting with exact length filtering completed successfully")
