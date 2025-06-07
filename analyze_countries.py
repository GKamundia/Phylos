#!/usr/bin/env python3
"""
Analyze countries in RVF metadata and update lat_longs.tsv
"""

import pandas as pd
import sys

def analyze_countries():
    """Analyze unique countries in the filtered metadata"""
    
    # Read the filtered metadata
    metadata_file = "results/filtered/rvf_metadata.tsv"
    try:
        df = pd.read_csv(metadata_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: Could not find {metadata_file}")
        return
    
    # Get unique countries
    countries = df['Country'].dropna().unique()
    countries_sorted = sorted(countries)
    
    print(f"Found {len(countries_sorted)} unique countries in the filtered metadata:")
    for i, country in enumerate(countries_sorted, 1):
        print(f"{i:2d}. {country}")
    
    # Read current lat_longs.tsv
    lat_longs_file = "config/lat_longs.tsv"
    try:
        lat_longs_df = pd.read_csv(lat_longs_file, sep='\t')
        existing_countries = set(lat_longs_df['location'].tolist())
    except FileNotFoundError:
        print(f"Error: Could not find {lat_longs_file}")
        return
    
    print(f"\nCountries already in lat_longs.tsv: {len(existing_countries)}")
    for country in sorted(existing_countries):
        print(f"  - {country}")
    
    # Find missing countries
    missing_countries = set(countries_sorted) - existing_countries
    
    print(f"\nMissing countries from lat_longs.tsv: {len(missing_countries)}")
    for country in sorted(missing_countries):
        print(f"  - {country}")
    
    return missing_countries, countries_sorted

if __name__ == "__main__":
    analyze_countries()
