#!/usr/bin/env python3
"""
Convert sequences.csv to proper FASTA and metadata TSV files for the Nextstrain pipeline
"""

import pandas as pd
import argparse
import os
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import sys
import re
from datetime import datetime

def setup_entrez(email, api_key=None):
    """Set up Entrez email and API key"""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

def fetch_sequence_by_accession(accession, retries=3):
    """Fetch sequence from NCBI by accession number"""
    for attempt in range(retries):
        try:
            print(f"Fetching sequence for {accession} (attempt {attempt + 1}/{retries})")
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return record
        except Exception as e:
            print(f"Failed to fetch {accession}: {e}")
            if attempt < retries - 1:
                time.sleep(2 ** attempt)  # Exponential backoff
            else:
                print(f"Failed to fetch {accession} after {retries} attempts")
                return None

def clean_strain_name(accession, isolate, organism):
    """Create a clean strain name from available information"""
    if isolate and isolate.strip():
        # Clean the isolate name
        strain = re.sub(r'[^\w\-_/.]', '_', isolate.strip())
    else:
        # Use accession as fallback
        strain = accession
    
    # Ensure strain doesn't start with a number or special character
    if strain and not strain[0].isalpha():
        strain = f"RVFV_{strain}"
    
    return strain

def clean_date(date_str):
    """Clean and standardize date format"""
    if not date_str or pd.isna(date_str):
        return ""
    
    date_str = str(date_str).strip()
    if not date_str:
        return ""
    
    # Try to parse various date formats
    date_patterns = [
        "%Y-%m-%d",
        "%Y/%m/%d", 
        "%Y",
        "%m/%d/%Y",
        "%d/%m/%Y"
    ]
    
    for pattern in date_patterns:
        try:
            parsed_date = datetime.strptime(date_str, pattern)
            return parsed_date.strftime("%Y-%m-%d")
        except ValueError:
            continue
    
    # If no pattern matches, try to extract just the year
    year_match = re.search(r'\b(19|20)\d{2}\b', date_str)
    if year_match:
        return year_match.group(0)
    
    return ""

def convert_csv_to_nextstrain_format(csv_file, output_fasta, output_metadata, email, api_key=None, segment_filter=None):
    """Convert sequences.csv to FASTA and metadata TSV files"""
    
    # Set up NCBI access
    setup_entrez(email, api_key)
    
    # Read the CSV file
    print(f"Reading {csv_file}...")
    df = pd.read_csv(csv_file)
    
    print(f"Found {len(df)} sequences in CSV")
    
    # Filter by segment if specified
    if segment_filter and segment_filter != "all":
        df = df[df['Segment'].str.upper() == segment_filter.upper()]
        print(f"Filtered to {len(df)} sequences for segment {segment_filter}")
    
    # Fetch sequences and create metadata
    sequences = []
    metadata_records = []
    
    for idx, row in df.iterrows():
        accession = row['Accession']
        print(f"Processing {idx+1}/{len(df)}: {accession}")
        
        # Fetch sequence
        seq_record = fetch_sequence_by_accession(accession)
        if not seq_record:
            print(f"Skipping {accession} - could not fetch sequence")
            continue
        
        # Create clean strain name
        strain = clean_strain_name(accession, row.get('Isolate', ''), row.get('Organism_Name', ''))
        
        # Update sequence record with clean name
        seq_record.id = strain
        seq_record.name = strain
        seq_record.description = f"{strain} | {accession}"
        
        sequences.append(seq_record)
        
        # Create metadata record
        metadata = {
            'strain': strain,
            'virus': 'Rift Valley fever virus',
            'accession': accession,
            'date': clean_date(row.get('Collection_Date', '')),
            'country': str(row.get('Country', '')).strip() if pd.notna(row.get('Country', '')) else '',
            'division': '',  # Not available in source data
            'location': str(row.get('Geo_Location', '')).strip() if pd.notna(row.get('Geo_Location', '')) else '',
            'host': str(row.get('Host', '')).strip() if pd.notna(row.get('Host', '')) else '',
            'segment': str(row.get('Segment', '')).strip() if pd.notna(row.get('Segment', '')) else '',
            'length': str(row.get('Length', '')).strip() if pd.notna(row.get('Length', '')) else str(len(seq_record.seq)),
            'Nuc_Completeness': str(row.get('Nuc_Completeness', '')).strip() if pd.notna(row.get('Nuc_Completeness', '')) else '',
            'latitude': '',  # Not available in source data
            'longitude': '',  # Not available in source data
        }
        
        metadata_records.append(metadata)
        
        # Add a small delay to be nice to NCBI
        time.sleep(0.1)
    
    # Write FASTA file
    print(f"Writing {len(sequences)} sequences to {output_fasta}")
    with open(output_fasta, 'w') as f:
        SeqIO.write(sequences, f, 'fasta')
    
    # Write metadata TSV file
    print(f"Writing metadata to {output_metadata}")
    metadata_df = pd.DataFrame(metadata_records)
    metadata_df.to_csv(output_metadata, sep='\t', index=False)
    
    print(f"Conversion complete!")
    print(f"- Sequences: {len(sequences)}")
    print(f"- Output FASTA: {output_fasta}")
    print(f"- Output metadata: {output_metadata}")
    
    return len(sequences)

def main():
    parser = argparse.ArgumentParser(description="Convert sequences.csv to Nextstrain format")
    parser.add_argument("--input-csv", required=True, help="Input sequences.csv file")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file")
    parser.add_argument("--output-metadata", required=True, help="Output metadata TSV file")
    parser.add_argument("--email", required=True, help="Email for NCBI access")
    parser.add_argument("--api-key", help="NCBI API key (optional)")
    parser.add_argument("--segment", choices=['L', 'M', 'S', 'all'], default='all', 
                        help="Filter by segment")
    
    args = parser.parse_args()
    
    # Create output directories
    os.makedirs(os.path.dirname(args.output_fasta), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_metadata), exist_ok=True)
    
    # Convert the data
    convert_csv_to_nextstrain_format(
        args.input_csv,
        args.output_fasta,
        args.output_metadata,
        args.email,
        args.api_key,
        args.segment
    )

if __name__ == "__main__":
    main()
