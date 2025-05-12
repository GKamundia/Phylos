#!/usr/bin/env python3
# filepath: rvf-nextstrain/scripts/prepare_metadata.py
"""
Clean and standardize metadata for Rift Valley Fever Nextstrain analysis
"""

import os
import pandas as pd
import argparse
import re
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(
        description="Clean and standardize RVF metadata")
    parser.add_argument('input', help="Path to input metadata TSV file")
    parser.add_argument('output', help="Path to output metadata TSV file")
    return parser.parse_args()

def standardize_dates(date_str):
    """Convert various date formats to YYYY-MM-DD"""
    if not date_str or pd.isna(date_str):
        return ""
    
    # Try common date formats
    date_formats = [
        "%Y-%m-%d", "%Y/%m/%d",  # YYYY-MM-DD, YYYY/MM/DD
        "%d-%m-%Y", "%d/%m/%Y",  # DD-MM-YYYY, DD/MM/YYYY
        "%m-%d-%Y", "%m/%d/%Y",  # MM-DD-YYYY, MM/DD/YYYY
        "%Y-%m", "%Y/%m",        # YYYY-MM, YYYY/MM
        "%Y"                      # YYYY
    ]
    
    # Try to parse the date
    for fmt in date_formats:
        try:
            date_obj = datetime.strptime(str(date_str), fmt)
            if fmt == "%Y":
                return f"{date_obj.year}"
            elif fmt in ["%Y-%m", "%Y/%m"]:
                return f"{date_obj.year}-{date_obj.month:02d}"
            else:
                return date_obj.strftime("%Y-%m-%d")
        except ValueError:
            continue
    
    # Extract year from string
    year_match = re.search(r'(\d{4})', str(date_str))
    if year_match:
        return year_match.group(1)
    
    return ""

def main():
    args = parse_args()
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Load metadata
    print(f"Loading metadata from {args.input}")
    try:
        metadata = pd.read_csv(args.input, sep='\t', dtype=str)
    except Exception as e:
        print(f"Error loading metadata: {e}")
        return
    
    print(f"Loaded {len(metadata)} records")
    
    # Clean and standardize
    print("Standardizing metadata...")
    
    # Standardize dates
    if 'date' in metadata.columns:
        metadata['date'] = metadata['date'].apply(standardize_dates)
        print(f"Standardized dates: {metadata['date'].nunique()} unique values")
    
    # Ensure strain names are unique
    if 'strain' in metadata.columns:
        # If duplicates exist, append a suffix
        duplicated = metadata['strain'].duplicated()
        if any(duplicated):
            print(f"Found {sum(duplicated)} duplicate strain names")
            dup_strains = metadata.loc[duplicated, 'strain'].unique()
            for strain in dup_strains:
                mask = metadata['strain'] == strain
                count = sum(mask)
                # Add suffix to duplicates
                indices = metadata.index[mask]
                for i, idx in enumerate(indices[1:], start=1):
                    metadata.loc[idx, 'strain'] = f"{strain}_{i}"
            print("Added suffix to duplicate strain names")
    
    # Ensure country is filled
    if 'country' in metadata.columns:
        empty_country = metadata['country'].isna() | (metadata['country'] == '')
        if any(empty_country):
            print(f"Warning: {sum(empty_country)} records with missing country information")
    
    # Save cleaned metadata
    print(f"Saving cleaned metadata to {args.output}")
    metadata.to_csv(args.output, sep='\t', index=False)
    print("Done!")

if __name__ == "__main__":
    main()