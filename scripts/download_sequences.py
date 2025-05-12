#!/usr/bin/env python3
"""
Download pathogen sequences from NCBI GenBank
This script is designed to be pathogen-agnostic, supporting various search terms
"""

import os
import yaml
import argparse
from Bio import Entrez
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="Download pathogen sequences from NCBI GenBank")
    parser.add_argument('--search-term', required=True,
                        help="NCBI search term for the pathogen")
    parser.add_argument('--max-sequences', type=int, default=500,
                        help="Maximum number of sequences to download")
    parser.add_argument('--output-sequences', required=True,
                        help="Path to output FASTA file")
    parser.add_argument('--output-metadata', required=True,
                        help="Path to output metadata TSV file")
    parser.add_argument('--email', required=True,
                        help="Email for NCBI Entrez (required by NCBI)")
    parser.add_argument('--segment', default=None,
                        help="Genome segment to filter for (if applicable)")
    return parser.parse_args()

def search_and_download(search_term, max_sequences, email, output_fasta, output_metadata, segment=None):
    """Search for and download sequences from NCBI"""
    print(f"Using email {email} for NCBI Entrez")
    Entrez.email = email
    
    # Create output directories if they don't exist
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    os.makedirs(os.path.dirname(output_metadata), exist_ok=True)
    
    print(f"Searching NCBI for: {search_term}")
    # Search for sequences
    search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=max_sequences)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    id_list = search_results["IdList"]
    print(f"Found {len(id_list)} sequences")
    
    if not id_list:
        print("No sequences found")
        return
    
    # Download sequences in batches
    batch_size = 100
    sequences = []
    metadata_lines = ["strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength"]
    
    for i in range(0, len(id_list), batch_size):
        batch_ids = id_list[i:i+batch_size]
        print(f"Downloading batch {i//batch_size + 1} of {(len(id_list)-1)//batch_size + 1}")
        fetch_handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="text")
        
        for record in SeqIO.parse(fetch_handle, "genbank"):
            # Process sequence
            sequences.append(record)
            
            # Extract metadata
            accession = record.id
            strain = record.annotations.get("strain", accession)
            length = len(record.seq)
            
            # Extract date
            date = ""
            for feature in record.features:
                if feature.type == "source":
                    if "collection_date" in feature.qualifiers:
                        date = feature.qualifiers["collection_date"][0]
                    break
            
            # Extract location
            country = ""
            division = ""
            for feature in record.features:
                if feature.type == "source":
                    if "country" in feature.qualifiers:
                        country = feature.qualifiers["country"][0]
                    break
            
            # Extract host
            host = ""
            for feature in record.features:
                if feature.type == "source":
                    if "host" in feature.qualifiers:
                        host = feature.qualifiers["host"][0]
                    break
            
            # Extract segment information
            segment_info = ""
            for feature in record.features:
                if feature.type == "source":
                    if "segment" in feature.qualifiers:
                        segment_info = feature.qualifiers["segment"][0]
                    break
            
            # Filter by segment if applicable
            if segment and segment_info != segment:
                continue
            
            # Add to metadata
            metadata_lines.append(f"{strain}\tPathogen\t{accession}\t{date}\t{country}\t{division}\t\t{host}\t{segment_info}\t{length}")
        
        fetch_handle.close()
    
    # Write sequences to FASTA file
    with open(output_fasta, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    
    # Write metadata to TSV file
    with open(output_metadata, "w") as f:
        f.write("\n".join(metadata_lines))
    
    print(f"Downloaded {len(sequences)} sequences")
    print(f"Sequences saved to: {output_fasta}")
    print(f"Metadata saved to: {output_metadata}")

def main():
    args = parse_args()
    
    search_term = args.search_term
    max_sequences = args.max_sequences
    email = args.email
    output_fasta = args.output_sequences
    output_metadata = args.output_metadata
    segment = args.segment
    
    search_and_download(search_term, max_sequences, email, output_fasta, output_metadata, segment)

if __name__ == "__main__":
    main()