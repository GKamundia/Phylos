#!/usr/bin/env python3
"""
Analyze segment data to understand reference sequence needs
"""

import pandas as pd
from Bio import SeqIO
import os

# Check filtered data for potential reference sequences
metadata_file = 'results/filtered/rvf_metadata.tsv'
sequences_file = 'results/filtered/rvf_filtered.fasta'

if os.path.exists(metadata_file):
    df = pd.read_csv(metadata_file, sep='\t')
    print('Segment distribution:')
    print(df['segment'].value_counts().sort_index())
    print('\nSample of data:')
    print(df[['strain', 'segment', 'Length', 'Country', 'date']].head(10))
    
    # Find sequences with complete dates for potential references
    complete_data = df[df['date'].notna() & (df['date'] != '')]
    print(f'\nSequences with complete dates: {len(complete_data)}')
    
    for segment in ['l', 'm', 's']:
        seg_data = complete_data[complete_data['segment'] == segment]
        if len(seg_data) > 0:
            # Find longest sequence for this segment as potential reference
            longest = seg_data.loc[seg_data['Length'].idxmax()]
            print(f'\n{segment.upper()} segment potential reference:')
            print(f'  Strain: {longest["strain"]}')
            print(f'  Length: {longest["Length"]} bp')
            print(f'  Date: {longest["date"]}')
            print(f'  Country: {longest["Country"]}')
            
    # Analyze length variations within segments
    print('\n=== SEGMENT LENGTH ANALYSIS ===')
    for segment in ['l', 'm', 's']:
        seg_data = df[df['segment'] == segment]
        if len(seg_data) > 0:
            lengths = seg_data['Length']
            print(f'\n{segment.upper()} segment:')
            print(f'  Count: {len(seg_data)}')
            print(f'  Length range: {lengths.min()} - {lengths.max()} bp')
            print(f'  Mean length: {lengths.mean():.1f} bp')
            print(f'  Length variation: {lengths.max() - lengths.min()} bp')
            print(f'  Coefficient of variation: {lengths.std()/lengths.mean()*100:.1f}%')

else:
    print('Filtered metadata not found')
