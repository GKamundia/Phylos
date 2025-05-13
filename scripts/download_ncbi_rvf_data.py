#!/usr/bin/env python3
"""
Script to download RVF virus sequences and metadata from NCBI
"""

import os
import argparse
import sys
from datetime import datetime
from Bio import Entrez, SeqIO
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Download RVF sequences from NCBI")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--api-key", help="NCBI API key for higher rate limits")
    parser.add_argument("--search-term", default="Rift Valley fever virus[orgn]", 
                        help="Search term for NCBI Entrez")
    parser.add_argument("--max-sequences", type=int, default=1000, 
                        help="Maximum number of sequences to download")
    parser.add_argument("--output-sequences", default="data/sequences/raw/rvf_sequences.fasta",
                        help="Output path for sequences FASTA file")
    parser.add_argument("--output-metadata", default="data/metadata/raw/rvf_metadata.tsv",
                        help="Output path for metadata TSV file")
    parser.add_argument("--segment", choices=['L', 'M', 'S', 'all'], default='all',
                        help="Filter by specific segment")
    return parser.parse_args()

def search_ncbi(search_term, max_results=1000, email=None):
    """Search NCBI for sequences matching the search term"""
    if email:
        Entrez.email = email
    else:
        Entrez.email = "example@example.com"
        
    print(f"Searching NCBI for: {search_term}")
    handle = Entrez.esearch(db="nucleotide", term=search_term, 
                           retmax=max_results, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record["IdList"]
    web_env = record["WebEnv"]
    query_key = record["QueryKey"]
    
    return id_list, web_env, query_key

def fetch_sequences(id_list, web_env=None, query_key=None, batch_size=100):
    """Fetch sequences by ID list"""
    sequences = []
    metadata_records = []
    segment_counts = {'L': 0, 'M': 0, 'S': 0, 'unknown': 0}
    batch_size = min(batch_size, len(id_list))
    
    print(f"Fetching {len(id_list)} sequences")
    for start in range(0, len(id_list), batch_size):
        end = min(start + batch_size, len(id_list))
        batch_ids = id_list[start:end]
        
        print(f"Fetching batch {start//batch_size + 1} ({start}-{end})")
        
        try:
            # Use history if provided
            if web_env and query_key:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    retstart=start,
                    retmax=batch_size,
                    webenv=web_env,
                    query_key=query_key
                )
            else:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    id=",".join(batch_ids)
                )
                
            # Parse GenBank records
            for record in SeqIO.parse(fetch_handle, "genbank"):
                accession = record.id
                
                # Extract metadata
                metadata = {
                    "strain": record.name or accession,
                    "virus": "Rift Valley fever virus",
                    "accession": accession,
                    "date": "",
                    "country": "",
                    "division": "",
                    "location": "",
                    "host": "",
                    "segment": "",
                    "length": len(record.seq),
                    "latitude": "",
                    "longitude": ""
                }
                
                # Extract details from features
                for feature in record.features:
                    if feature.type == "source":
                        # Date
                        if "collection_date" in feature.qualifiers:
                            metadata["date"] = feature.qualifiers["collection_date"][0]
                        
                        # Country
                        if "country" in feature.qualifiers:
                            country_field = feature.qualifiers["country"][0]
                            # Handle potential division information (country: division)
                            if ":" in country_field:
                                country, division = country_field.split(":", 1)
                                metadata["country"] = country.strip()
                                metadata["division"] = division.strip()
                            else:
                                metadata["country"] = country_field.strip()
                        
                        # Host
                        if "host" in feature.qualifiers:
                            metadata["host"] = feature.qualifiers["host"][0]
                        
                        # Segment identification
                        segment_found = False
                        for qualifier in ["product", "note", "segment"]:
                            if qualifier in feature.qualifiers:
                                value = feature.qualifiers[qualifier][0].lower()
                                if "segment l" in value or "large" in value or "rdrp" in value:
                                    metadata["segment"] = "L"
                                    segment_counts["L"] += 1
                                    segment_found = True
                                    break
                                elif "segment m" in value or "medium" in value or "glyco" in value:
                                    metadata["segment"] = "M"
                                    segment_counts["M"] += 1
                                    segment_found = True
                                    break
                                elif "segment s" in value or "small" in value or "nucleo" in value:
                                    metadata["segment"] = "S"
                                    segment_counts["S"] += 1
                                    segment_found = True
                                    break
                
                # Default to unknown if segment not identified
                if not segment_found:
                    metadata["segment"] = "unknown"
                    segment_counts["unknown"] += 1
                
                # Add to our collections
                sequences.append(record)
                metadata_records.append(metadata)
                
            fetch_handle.close()
            
        except Exception as e:
            print(f"Error fetching batch {start}-{end}: {e}")
            continue
    
    print(f"Segment counts: L={segment_counts['L']}, M={segment_counts['M']}, "
          f"S={segment_counts['S']}, unknown={segment_counts['unknown']}")
    
    return sequences, metadata_records

def main():
    args = parse_args()
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(args.output_sequences), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_metadata), exist_ok=True)
    
    # Set up Entrez API key if provided
    if args.api_key:
        Entrez.api_key = args.api_key
    
    # Search NCBI
    id_list, web_env, query_key = search_ncbi(
        args.search_term, args.max_sequences, args.email
    )
    
    if not id_list:
        print("No sequences found matching the search term")
        return 1
    
    # Fetch sequences and metadata
    sequences, metadata_records = fetch_sequences(id_list, web_env, query_key)
    
    if not sequences:
        print("Failed to fetch any sequences")
        return 1
    
    # Filter by segment if requested
    if args.segment != 'all':
        filtered_sequences = []
        filtered_metadata = []
        
        for seq, meta in zip(sequences, metadata_records):
            if meta["segment"] == args.segment:
                filtered_sequences.append(seq)
                filtered_metadata.append(meta)
        
        sequences = filtered_sequences
        metadata_records = filtered_metadata
        print(f"Filtered to {len(sequences)} sequences of segment {args.segment}")
    
    # Write sequences to FASTA
    with open(args.output_sequences, 'w') as f:
        SeqIO.write(sequences, f, 'fasta')
    
    # Write metadata to TSV
    if metadata_records:
        df = pd.DataFrame(metadata_records)
        df.to_csv(args.output_metadata, sep='\t', index=False)
    else:
        # Create empty metadata file with headers
        with open(args.output_metadata, 'w') as f:
            headers = ["strain", "virus", "accession", "date", "country", 
                      "division", "location", "host", "segment", "length",
                      "latitude", "longitude"]
            f.write('\t'.join(headers) + '\n')
    
    print(f"Downloaded {len(sequences)} sequences to {args.output_sequences}")
    print(f"Wrote metadata to {args.output_metadata}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())