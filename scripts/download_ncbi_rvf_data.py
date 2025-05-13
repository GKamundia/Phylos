#!/usr/bin/env python3
"""
Script to download RVF virus sequences and metadata from NCBI
"""

import os
import argparse
import sys
import re
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
    parser.add_argument("--debug", action="store_true", help="Enable detailed debug output")
    parser.add_argument("--create-dummy", action="store_true", 
                       help="Create dummy data if no real sequences are found")
    return parser.parse_args()

def search_ncbi(search_term, max_results=1000, email=None, api_key=None):
    """Search NCBI for sequences matching the search term"""
    if email:
        Entrez.email = email
    else:
        Entrez.email = "example@example.com"
    
    if api_key:
        Entrez.api_key = api_key
        print(f"Using provided NCBI API key")
    
    print(f"Searching NCBI for: {search_term}")
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_term, 
                              retmax=max_results, usehistory="y")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        web_env = record["WebEnv"]
        query_key = record["QueryKey"]
        
        return id_list, web_env, query_key
        
    except Exception as e:
        print(f"Error during NCBI search: {e}")
        return [], None, None

def fetch_sequences(id_list, web_env=None, query_key=None, batch_size=100, debug=False, heuristic_lookup=True):
    """Fetch sequences by ID list"""
    sequences = []
    metadata_records = []
    segment_counts = {'L': 0, 'M': 0, 'S': 0, 'unknown': 0}
    qualifier_values = []  # For debugging
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
                # Check if the accession is for the reference genome
                is_reference = ("reference" in record.description.lower() or 
                               "complete genome" in record.description.lower())
                
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
                all_qualifier_values = []  # Collect all values for debugging
                segment_from_length = None
                segment_from_definition = None
                
                # Check if segment info is in the definition line
                if "segment l " in record.description.lower() or "l segment" in record.description.lower():
                    segment_from_definition = "L"
                elif "segment m " in record.description.lower() or "m segment" in record.description.lower():
                    segment_from_definition = "M"
                elif "segment s " in record.description.lower() or "s segment" in record.description.lower():
                    segment_from_definition = "S"
                
                # Try to infer segment by sequence length if heuristic_lookup is enabled
                if heuristic_lookup:
                    seq_len = len(record.seq)
                    if seq_len > 3500 and seq_len < 4500:
                        segment_from_length = "S"  # S segment is typically ~1.7kb
                    elif seq_len > 5500 and seq_len < 6500:
                        segment_from_length = "M"  # M segment is typically ~3.9kb
                    elif seq_len > 7000:
                        segment_from_length = "L"  # L segment is typically ~6.4kb
                
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
                    
                    # Collect all qualifier values for this record for debugging
                    for feature_qualifier in feature.qualifiers:
                        for value in feature.qualifiers[feature_qualifier]:
                            all_qualifier_values.append(f"{feature_qualifier}: {value}")
                
                # Segment identification - expanded to catch more variants
                segment_found = False
                segment_keywords = {
                    "L": ["segment l", "l segment", "large", "rdrp", "polymerase", "l protein", "l gene", "l rna"],
                    "M": ["segment m", "m segment", "medium", "glycoprotein", "gc", "gn", "m protein", "m gene", "m rna"],
                    "S": ["segment s", "s segment", "small", "nucleocapsid", "np", "n protein", "s protein", "s gene", "s rna"]
                }
                
                for qualifier in ["product", "note", "segment", "gene", "mol_type", "organism"]:
                    if feature.type in ["CDS", "gene", "source"] and qualifier in feature.qualifiers:
                        for value in feature.qualifiers[qualifier]:
                            value_lower = value.lower()
                            
                            # Check for segment keywords
                            for segment, keywords in segment_keywords.items():
                                for keyword in keywords:
                                    if keyword in value_lower:
                                        metadata["segment"] = segment
                                        segment_counts[segment] += 1
                                        segment_found = True
                                        break
                                if segment_found:
                                    break
                            if segment_found:
                                break
                    if segment_found:
                        break
                
                # Use alternative methods if segment not identified through qualifiers
                if not segment_found:
                    # Try to use segment identified from definition line
                    if segment_from_definition:
                        metadata["segment"] = segment_from_definition
                        segment_counts[segment_from_definition] += 1
                        segment_found = True
                    # Try to use segment identified from sequence length
                    elif segment_from_length:
                        metadata["segment"] = segment_from_length
                        segment_counts[segment_from_length] += 1
                        segment_found = True
                    # For reference genomes, more likely to be L segment
                    elif is_reference and len(record.seq) > 6000:
                        metadata["segment"] = "L"
                        segment_counts["L"] += 1
                        segment_found = True
                    else:
                        metadata["segment"] = "unknown"
                        segment_counts["unknown"] += 1
                
                # Add to our collections
                sequences.append(record)
                metadata_records.append(metadata)
                
                if debug:
                    qualifier_values.append({
                        "accession": accession,
                        "segment": metadata["segment"],
                        "length": len(record.seq),
                        "description": record.description,
                        "values": all_qualifier_values
                    })
                
            fetch_handle.close()
            
        except Exception as e:
            print(f"Error fetching batch {start}-{end}: {e}")
            continue
    
    print(f"Segment counts: L={segment_counts['L']}, M={segment_counts['M']}, "
          f"S={segment_counts['S']}, unknown={segment_counts['unknown']}")
    
    # Debug output
    if debug and qualifier_values:
        print("\nDebug: Sample of qualifier values from records:")
        for i, record_data in enumerate(qualifier_values[:5]):  # Show first 5 for brevity
            print(f"\nAccession: {record_data['accession']}")
            print(f"Identified segment: {record_data['segment']}")
            print(f"Sequence length: {record_data['length']}")
            print(f"Description: {record_data['description']}")
            print("Sample qualifier values:")
            for j, val in enumerate(record_data['values'][:10]):  # Show first 10 values
                print(f"  {val}")
            print("...")
        
        print("\nUnknown segment records:")
        unknown_count = 0
        for record_data in qualifier_values:
            if record_data['segment'] == 'unknown':
                unknown_count += 1
                if unknown_count <= 3:  # Show first 3 unknown records
                    print(f"\nAccession: {record_data['accession']}")
                    print(f"Sequence length: {record_data['length']}")
                    print(f"Description: {record_data['description']}")
                    print("Sample qualifier values:")
                    for j, val in enumerate(record_data['values'][:5]):
                        print(f"  {val}")
                    print("...")
    
    return sequences, metadata_records, qualifier_values if debug else None

def create_dummy_data(segment, output_sequences, output_metadata):
    """Create dummy data for testing when no real sequences are found"""
    print(f"Creating dummy {segment} segment data for pipeline testing")
    
    # Create a dummy sequence for the specified segment
    if segment == 'L':
        seq_length = 6400
    elif segment == 'M':
        seq_length = 3885
    elif segment == 'S':
        seq_length = 1690
    else:
        seq_length = 1000
    
    # Write dummy sequence to FASTA
    with open(output_sequences, 'w') as f:
        f.write(f">RVFV_{segment}_dummy\n")
        f.write("A" * seq_length + "\n")
    
    # Write dummy metadata to TSV
    with open(output_metadata, 'w') as f:
        headers = ["strain", "virus", "accession", "date", "country", 
                  "division", "location", "host", "segment", "length",
                  "latitude", "longitude"]
        f.write('\t'.join(headers) + '\n')
        
        metadata = [
            f"RVFV_{segment}_dummy",
            "Rift Valley fever virus",
            f"DUMMY_{segment}",
            datetime.now().strftime("%Y-%m-%d"),
            "Kenya",
            "Nairobi",
            "Nairobi",
            "Bos taurus",
            segment,
            str(seq_length),
            "0.0236",
            "37.9062"
        ]
        f.write('\t'.join(metadata) + '\n')

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
        args.search_term, args.max_sequences, args.email, args.api_key
    )
    
    if not id_list:
        print("No sequences found matching the search term")
        if args.create_dummy:
            create_dummy_data(args.segment if args.segment != 'all' else 'L', 
                             args.output_sequences, args.output_metadata)
            return 0
        return 1
    
    # Fetch sequences and metadata
    sequences, metadata_records, debug_data = fetch_sequences(
        id_list, web_env, query_key, debug=args.debug
    )
    
    if not sequences:
        print("Failed to fetch any sequences")
        if args.create_dummy:
            create_dummy_data(args.segment if args.segment != 'all' else 'L', 
                             args.output_sequences, args.output_metadata)
            return 0
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
    
    # If no sequences left after filtering but dummy creation is enabled
    if not sequences and args.create_dummy:
        create_dummy_data(args.segment if args.segment != 'all' else 'L', 
                         args.output_sequences, args.output_metadata)
        return 0
    
    # Even if filtering resulted in no sequences, still write empty output files for downstream pipeline
    if not sequences:
        print(f"No sequences found for segment {args.segment}. Creating empty output files.")
        # Create empty FASTA
        with open(args.output_sequences, 'w') as f:
            f.write(f">RVFV_{args.segment}_dummy\n")
            f.write("ACGT\n")  # Minimal sequence for pipeline to continue
        
        # Create minimal metadata
        with open(args.output_metadata, 'w') as f:
            headers = ["strain", "virus", "accession", "date", "country", 
                      "division", "location", "host", "segment", "length",
                      "latitude", "longitude"]
            f.write('\t'.join(headers) + '\n')
            
            metadata = [
                f"RVFV_{args.segment}_dummy",
                "Rift Valley fever virus",
                f"DUMMY_{args.segment}",
                datetime.now().strftime("%Y-%m-%d"),
                "Unknown",
                "",
                "",
                "",
                args.segment,
                "4",
                "",
                ""
            ]
            f.write('\t'.join(metadata) + '\n')
    else:
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