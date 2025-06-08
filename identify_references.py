#!/usr/bin/env python3
"""
Identify potential reference sequences for each RVF segment
"""

import pandas as pd
from Bio import SeqIO
import os

def main():
    # Load filtered metadata
    metadata_path = "results/filtered/rvf_metadata.tsv"
    sequences_path = "results/filtered/rvf_filtered.fasta"
    
    if not os.path.exists(metadata_path) or not os.path.exists(sequences_path):
        print("Filtered data not found. Run filtering first.")
        return
    
    df = pd.read_csv(metadata_path, sep='\t')
    
    # Filter for complete sequences only
    complete = df[df['Nuc_Completeness'] == 'complete']
    print(f"Total complete sequences: {len(complete)}")
    print(f"Segment distribution:")
    print(complete['segment'].value_counts())
    
    # Load sequences
    sequences = {record.id: record for record in SeqIO.parse(sequences_path, "fasta")}
    
    # Find best reference candidates for each segment
    segment_refs = {}
    
    for segment in ['l', 'm', 's']:
        seg_data = complete[complete['segment'] == segment]
        
        if len(seg_data) == 0:
            print(f"\nNo complete sequences for segment {segment.upper()}")
            continue
            
        print(f"\n{segment.upper()} Segment Analysis:")
        print(f"  Complete sequences: {len(seg_data)}")
        
        # Sort by length (descending) and date (recent first)
        seg_data_sorted = seg_data.sort_values(['Length', 'date'], ascending=[False, False])
        
        # Get top candidates
        top_candidates = seg_data_sorted.head(5)
        
        print(f"  Top {len(top_candidates)} candidates:")
        for idx, (_, row) in enumerate(top_candidates.iterrows()):
            strain = row['strain']
            seq_id = row['accession']
            
            # Check if sequence exists in FASTA
            if seq_id in sequences:
                actual_length = len(sequences[seq_id].seq)
                print(f"    {idx+1}. {strain} ({seq_id})")
                print(f"       Length: {row['Length']} bp (actual: {actual_length} bp)")
                print(f"       Date: {row['date']}")
                print(f"       Country: {row['country']}")
                
                # Select first available as reference
                if idx == 0:
                    segment_refs[segment] = {
                        'strain': strain,
                        'accession': seq_id,
                        'length': actual_length,
                        'date': row['date'],
                        'country': row['country']
                    }
            else:
                print(f"    {idx+1}. {strain} ({seq_id}) - SEQUENCE NOT FOUND")
    
    print(f"\n=== SELECTED REFERENCES ===")
    for segment, ref_info in segment_refs.items():
        print(f"{segment.upper()} segment: {ref_info['strain']} ({ref_info['accession']})")
        print(f"  Length: {ref_info['length']} bp, Date: {ref_info['date']}, Country: {ref_info['country']}")
    
    return segment_refs

if __name__ == "__main__":
    refs = main()
