#!/usr/bin/env python3
"""
Split metadata by segment for RVF Nextstrain pipeline
Creates segment-specific metadata files needed for export rules
"""

import pandas as pd
import argparse
import os
import sys

def split_metadata_by_segment(input_file, output_dir):
    """
    Split metadata file by segment and create segment-specific files
    
    Args:
        input_file: Path to input metadata TSV file
        output_dir: Directory to output segment-specific files
    
    Returns:
        Dict with segment names and number of records written
    """
    # Load metadata
    print(f"Loading metadata from {input_file}")
    metadata = pd.read_csv(input_file, sep='\t')
    print(f"Loaded {len(metadata)} records")
    
    # Check segment distribution
    if 'segment' not in metadata.columns:
        print("ERROR: No 'segment' column found in metadata")
        return {}
    
    segment_dist = metadata['segment'].value_counts()
    print(f"Segment distribution: {segment_dist.to_dict()}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Split by segment
    results = {}
    for segment in ['l', 'm', 's']:
        segment_data = metadata[metadata['segment'] == segment].copy()
        
        if len(segment_data) > 0:
            # Create output filename
            base_filename = os.path.basename(input_file)
            if base_filename.endswith('_metadata.tsv'):
                prefix = base_filename.replace('_metadata.tsv', '')
                output_file = os.path.join(output_dir, f"{prefix}_metadata_{segment.upper()}.tsv")
            else:
                output_file = os.path.join(output_dir, f"metadata_{segment.upper()}.tsv")
            
            # Write segment-specific metadata
            segment_data.to_csv(output_file, sep='\t', index=False)
            results[segment.upper()] = len(segment_data)
            print(f"Wrote {len(segment_data)} records for segment {segment.upper()} to {output_file}")
        else:
            print(f"No data for segment {segment.upper()}")
            results[segment.upper()] = 0
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Split metadata by segment")
    parser.add_argument("--input", required=True, help="Input metadata TSV file")
    parser.add_argument("--output-dir", required=True, help="Output directory for segment-specific files")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"ERROR: Input file {args.input} not found")
        return 1
    
    try:
        results = split_metadata_by_segment(args.input, args.output_dir)
        
        total_written = sum(results.values())
        print(f"\nSummary:")
        print(f"Total records written: {total_written}")
        for segment, count in results.items():
            print(f"  {segment}: {count} records")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
