#!/usr/bin/env python3
"""
Run filtering for RVF virus sequences to keep only complete sequences
"""

import sys
import os
import argparse
import pandas as pd
from Bio import SeqIO

def main():
    print("Starting run_filter.py script")
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--metadata", required=True) 
    parser.add_argument("--output-sequences", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()
    
    print(f"Arguments parsed: sequences={args.sequences}, metadata={args.metadata}, output={args.output_sequences}")
    
    # Create output directory
    os.makedirs(os.path.dirname(args.output_sequences), exist_ok=True)
    print(f"Created output directory: {os.path.dirname(args.output_sequences)}")
    
    log_content = []
    log_content.append("Filtering for complete sequences across all segments\n")
    
    try:
        # Read metadata
        print(f"Reading metadata from {args.metadata}")
        metadata = pd.read_csv(args.metadata, sep='\t')
        total_count = len(metadata)
        log_content.append(f"Total sequences before filtering: {total_count}\n")
        print(f"Total sequences before filtering: {total_count}")
        
        # Check if 'Nuc_Completeness' column exists
        if 'Nuc_Completeness' not in metadata.columns:
            error_msg = f"Error: 'Nuc_Completeness' column not found in metadata. Available columns: {list(metadata.columns)}"
            print(error_msg)
            log_content.append(f"{error_msg}\n")
            raise ValueError(error_msg)
          # Filter metadata to keep only complete sequences
        print("Filtering metadata for complete sequences")
        filtered_metadata = metadata[metadata['Nuc_Completeness'] == 'complete']
        kept_count_complete = len(filtered_metadata)
        dropped_count_complete = total_count - kept_count_complete
        
        log_content.append(f"Sequences after completeness filtering: {kept_count_complete}\n")
        log_content.append(f"Sequences dropped due to 'Nuc_Completeness != complete': {dropped_count_complete}\n")
        print(f"After completeness filter - retained: {kept_count_complete}, dropped: {dropped_count_complete}")
        
        # Additional filtering for missing country and date
        print("Applying additional filters for missing country and date")
        initial_after_complete = len(filtered_metadata)
        
        # Filter out records with missing country
        if 'country' in filtered_metadata.columns:
            before_country = len(filtered_metadata)
            filtered_metadata = filtered_metadata[
                (filtered_metadata['country'].notna()) & 
                (filtered_metadata['country'] != '') & 
                (filtered_metadata['country'] != 'nan')
            ]
            after_country = len(filtered_metadata)
            dropped_country = before_country - after_country
            log_content.append(f"Sequences dropped due to missing country: {dropped_country}\n")
            print(f"Country filter - retained: {after_country}, dropped: {dropped_country}")
        else:
            log_content.append("Warning: 'country' column not found, skipping country filter\n")
            print("Warning: 'country' column not found, skipping country filter")
        
        # Filter out records with missing date
        if 'date' in filtered_metadata.columns:
            before_date = len(filtered_metadata)
            filtered_metadata = filtered_metadata[
                (filtered_metadata['date'].notna()) & 
                (filtered_metadata['date'] != '') & 
                (filtered_metadata['date'] != 'nan')
            ]
            after_date = len(filtered_metadata)
            dropped_date = before_date - after_date
            log_content.append(f"Sequences dropped due to missing date: {dropped_date}\n")
            print(f"Date filter - retained: {after_date}, dropped: {dropped_date}")
        else:
            log_content.append("Warning: 'date' column not found, skipping date filter\n")
            print("Warning: 'date' column not found, skipping date filter")
        
        kept_count = len(filtered_metadata)
        total_dropped = total_count - kept_count
        log_content.append(f"Final sequences retained after all filtering: {kept_count}\n")
        log_content.append(f"Total sequences dropped: {total_dropped}\n")
        print(f"Final results - retained: {kept_count}, total dropped: {total_dropped}")
          # Get the list of sequence IDs to keep (use Accession which matches FASTA headers)
        if 'Accession' in filtered_metadata.columns:
            keep_ids = set(filtered_metadata['Accession'])
            print(f"Keeping {len(keep_ids)} unique accession IDs")
        else:
            # Fallback to strain if Accession not available
            keep_ids = set(filtered_metadata['strain'])
            print(f"Keeping {len(keep_ids)} unique strain IDs (fallback)")
        
        # Write filtered metadata
        print(f"Writing filtered metadata to {args.output_metadata}")
        filtered_metadata.to_csv(args.output_metadata, sep='\t', index=False)
        
        # Filter sequences
        sequence_count = 0
        print(f"Reading sequences from {args.sequences}")
        print(f"Writing filtered sequences to {args.output_sequences}")
        with open(args.sequences, "r") as input_handle, open(args.output_sequences, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                if record.id in keep_ids:
                    SeqIO.write(record, output_handle, "fasta")
                    sequence_count += 1
        
        log_content.append(f"Filtering complete. Written {sequence_count} sequences to {args.output_sequences}\n")
        print(f"Filtering complete. Written {sequence_count} sequences to {args.output_sequences}")
        
        # Segment distribution
        if 'Segment' in filtered_metadata.columns:
            segment_counts = filtered_metadata['Segment'].value_counts().to_dict()
            log_content.append("Segment distribution in filtered data:\n")
            print("Segment distribution in filtered data:")
            for segment, count in segment_counts.items():
                log_content.append(f"  {segment}: {count}\n")
                print(f"  {segment}: {count}")
        
        # Write log content to file
        print(f"Writing log to {args.log}")
        with open(args.log, 'w') as log_file:
            log_file.writelines(log_content)
        
    except Exception as e:
        print(f"Error during filtering: {e}", file=sys.stderr)
        # Try to write the error to the log file
        try:
            with open(args.log, 'w') as log_file:
                log_file.writelines(log_content)
                log_file.write(f"Error during filtering: {e}\n")
        except:
            print(f"Also failed to write to log file: {args.log}", file=sys.stderr)
        return 1
    
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

if __name__ == "__main__":
    sys.exit(main())
