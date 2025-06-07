#!/usr/bin/env python3
"""
Post-processing script to fix missing metadata fields in RVF data
This script adds the missing fields (Nuc_Completeness, Segment, etc.) 
based on sequence information and GenBank titles.
"""

import pandas as pd
import re
import argparse
from pathlib import Path

def extract_segment_from_title(title):
    """Extract segment from GenBank title"""
    title_lower = title.lower()
    
    if 'segment s' in title_lower:
        return 'S'
    elif 'segment m' in title_lower:
        return 'M'
    elif 'segment l' in title_lower:
        return 'L'
    
    return ''

def determine_completeness_from_title_and_length(title, length, segment):
    """Determine completeness from title and sequence length"""
    title_lower = title.lower()
    
    # Check title first
    if 'complete genome' in title_lower or 'complete sequence' in title_lower:
        return 'complete'
    elif 'partial' in title_lower:
        return 'partial'
    
    # Use length-based heuristic for known segments
    if segment and length > 0:
        # Known complete lengths for RVF segments (approximate)
        complete_lengths = {'S': 1690, 'M': 3885, 'L': 6404}
        
        if segment in complete_lengths:
            expected = complete_lengths[segment]
            # Consider complete if within 15% of expected length
            if abs(length - expected) / expected < 0.15:
                return 'complete'
    
    return 'partial'

def extract_geo_location_from_title(title):
    """Extract geographic location from title if possible"""
    # Look for common geographic patterns in titles
    geo_patterns = [
        r'(?:from|isolate)\s+([A-Z][a-z]+(?:\s+[A-Z][a-z]+)?)',  # Country names
        r'([A-Z][a-z]+)\s*:\s*[A-Z]',  # Country: region pattern
    ]
    
    for pattern in geo_patterns:
        match = re.search(pattern, title)
        if match:
            return match.group(1)
    
    return ''

def fix_genbank_refseq(accession):
    """Determine if accession is GenBank or RefSeq"""
    if accession.startswith(('NC_', 'NM_', 'NP_', 'XM_', 'XP_', 'XR_', 'NR_')):
        return 'RefSeq'
    else:
        return 'GenBank'

def extract_country_from_geo(geo_location):
    """Extract country from geo_location field"""
    if not geo_location:
        return ''
    
    if ':' in geo_location:
        return geo_location.split(':')[0].strip()
    else:
        return geo_location.strip()

def post_process_metadata(input_file, output_file=None):
    """Post-process metadata to add missing fields"""
    
    print(f"Reading metadata from: {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"Processing {len(df)} records...")
    
    # Fix each record
    for idx, row in df.iterrows():
        title = str(row.get('GenBank_Title', ''))
        length = int(row.get('Length', 0))
        accession = str(row.get('Accession', ''))
        
        # Extract/fix segment
        if pd.isna(row.get('Segment')) or not row.get('Segment'):
            segment = extract_segment_from_title(title)
            df.at[idx, 'Segment'] = segment
        else:
            segment = row.get('Segment')
        
        # Extract/fix completeness
        if pd.isna(row.get('Nuc_Completeness')) or not row.get('Nuc_Completeness'):
            completeness = determine_completeness_from_title_and_length(title, length, segment)
            df.at[idx, 'Nuc_Completeness'] = completeness
        
        # Fix GenBank/RefSeq
        if pd.isna(row.get('GenBank_RefSeq')) or not row.get('GenBank_RefSeq'):
            source = fix_genbank_refseq(accession)
            df.at[idx, 'GenBank_RefSeq'] = source
        
        # Extract geo location if missing
        if pd.isna(row.get('Geo_Location')) or not row.get('Geo_Location'):
            geo = extract_geo_location_from_title(title)
            df.at[idx, 'Geo_Location'] = geo
        
        # Extract country from geo_location
        geo_location = df.at[idx, 'Geo_Location']
        if geo_location and (pd.isna(row.get('Country')) or not row.get('Country')):
            country = extract_country_from_geo(geo_location)
            df.at[idx, 'Country'] = country
    
    # Print summary
    print(f"\n=== POST-PROCESSING SUMMARY ===")
    print(f"Total records: {len(df)}")
    
    if 'Segment' in df.columns:
        segment_counts = df['Segment'].value_counts()
        print(f"Segments: {dict(segment_counts)}")
    
    if 'Nuc_Completeness' in df.columns:
        completeness_counts = df['Nuc_Completeness'].value_counts()
        print(f"Completeness: {dict(completeness_counts)}")
    
    if 'GenBank_RefSeq' in df.columns:
        source_counts = df['GenBank_RefSeq'].value_counts()
        print(f"Sources: {dict(source_counts)}")
    
    # Save result
    if output_file is None:
        output_file = input_file.replace('.tsv', '_fixed.tsv')
    
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Fixed metadata saved to: {output_file}")
    
    return df

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Post-process RVF metadata to fix missing fields")
    parser.add_argument("input_file", help="Input metadata TSV file")
    parser.add_argument("--output", help="Output file path (default: input_file_fixed.tsv)")
    
    args = parser.parse_args()
    
    if not Path(args.input_file).exists():
        print(f"Error: Input file {args.input_file} not found")
        return
    
    post_process_metadata(args.input_file, args.output)

if __name__ == "__main__":
    main()
