#!/usr/bin/env python3
import pandas as pd
import json
import os

print('=== STAGE 2: DATA FILTERING AND METADATA PREPARATION COMPLETE ===')
print()

# Read prepared metadata
prepared_df = pd.read_csv('data/metadata/rvf_metadata.tsv', sep='\t')
print(f'Prepared metadata records: {len(prepared_df)}')
print(f'Prepared metadata columns: {len(prepared_df.columns)}')

# Read filtered data
filtered_df = pd.read_csv('results/filtered/rvf_metadata.tsv', sep='\t')
print(f'Filtered metadata records: {len(filtered_df)}')

# Read validation report
with open('data/metadata/rvf_validation_report.json', 'r') as f:
    validation = json.load(f)

print()
print('=== METADATA PREPARATION RESULTS ===')
print(f'Valid records: {validation["valid"]}')
print(f'Total records processed: {validation["total_records"]}')
print(f'Invalid records: {validation["invalid_records"]}')

print()
print('=== FILTERING RESULTS ===')
print(f'Input sequences: 1,811')
print(f'Filtered sequences: {len(filtered_df)} ({len(filtered_df)/1811*100:.1f}%)')
print(f'Sequences dropped: {1811-len(filtered_df)} ({(1811-len(filtered_df))/1811*100:.1f}%)')

print()
print('=== SEGMENT DISTRIBUTION IN FILTERED DATA ===')
if 'Segment' in filtered_df.columns:
    segment_counts = filtered_df['Segment'].value_counts()
    for segment, count in segment_counts.items():
        print(f'  {segment}: {count} sequences')

print()
print('=== GEOGRAPHIC DISTRIBUTION ===')
if 'country' in filtered_df.columns:
    country_counts = filtered_df['country'].value_counts().head(10)
    print('Top 10 countries:')
    for country, count in country_counts.items():
        print(f'  {country}: {count} sequences')

print()
print('=== DATE COVERAGE ===')
if 'date' in filtered_df.columns:
    date_stats = filtered_df['date'].str.len().value_counts().sort_index()
    print('Date completeness:')
    for length, count in date_stats.items():
        if length == 10:
            print(f'  Complete dates (YYYY-MM-DD): {count}')
        elif length == 4:
            print(f'  Year-only dates (YYYY): {count}')
        else:
            print(f'  Other date formats (length {length}): {count}')

print()
print('=== COORDINATE COVERAGE ===')
if 'latitude' in filtered_df.columns and 'longitude' in filtered_df.columns:
    coords_available = ((filtered_df['latitude'].notna()) & (filtered_df['longitude'].notna())).sum()
    print(f'Records with coordinates: {coords_available} ({coords_available/len(filtered_df)*100:.1f}%)')
    print(f'Records missing coordinates: {len(filtered_df)-coords_available} ({(len(filtered_df)-coords_available)/len(filtered_df)*100:.1f}%)')

print()
print('=== FILES CREATED ===')
print(f'✅ Raw sequences: data/sequences/raw/rvf_sequences.fasta ({os.path.getsize("data/sequences/raw/rvf_sequences.fasta") / (1024*1024):.1f} MB)')
print(f'✅ Raw metadata: data/metadata/raw/rvf_metadata.tsv ({os.path.getsize("data/metadata/raw/rvf_metadata.tsv") / 1024:.1f} KB)')
print(f'✅ Prepared metadata: data/metadata/rvf_metadata.tsv ({os.path.getsize("data/metadata/rvf_metadata.tsv") / 1024:.1f} KB)')
print(f'✅ Filtered sequences: results/filtered/rvf_filtered.fasta ({os.path.getsize("results/filtered/rvf_filtered.fasta") / (1024*1024):.1f} MB)')
print(f'✅ Filtered metadata: results/filtered/rvf_metadata.tsv ({os.path.getsize("results/filtered/rvf_metadata.tsv") / 1024:.1f} KB)')
print(f'✅ Validation report: data/metadata/rvf_validation_report.json')
