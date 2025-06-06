#!/usr/bin/env python3
"""
Script to download RVF virus sequences and metadata from NCBI Virus
This script mimics the exact search and data retrieval used by the NCBI Virus website
to ensure perfect data matching including all segments and correct species names.
"""

import os
import sys
import json
import time
import re
from datetime import datetime
from Bio import Entrez, SeqIO
import pandas as pd

def setup_entrez(email, api_key=None):
    """Setup Entrez with email and optional API key"""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        print(f"Using NCBI API key for enhanced access")

def search_ncbi_virus_exact(email, api_key=None, max_sequences=None):
    """
    Search NCBI using the exact same approach as the NCBI Virus website
    This ensures we get ALL RVF virus sequences including all segments
    """
    setup_entrez(email, api_key)
    
    # Multiple search strategies to ensure we get ALL RVF sequences
    search_strategies = [
        # Primary search - matches NCBI Virus website exactly
        "txid11588[Organism:exp]",  # Exact taxid search like the website
        "Phlebovirus riftense[Organism]",  # Updated species name
        "Rift Valley fever virus[Organism]",  # Original name (for older sequences)
        # Additional searches to catch all variants
        "Rift Valley fever virus[Title] AND (segment OR complete OR partial)",
        "RVFV[Title] AND virus[Title]",
        # Segment-specific searches to ensure we get all segments
        "(Rift Valley fever virus OR Phlebovirus riftense) AND segment L",
        "(Rift Valley fever virus OR Phlebovirus riftense) AND segment M", 
        "(Rift Valley fever virus OR Phlebovirus riftense) AND segment S",
        "(Rift Valley fever virus OR Phlebovirus riftense) AND (large segment OR L segment)",
        "(Rift Valley fever virus OR Phlebovirus riftense) AND (medium segment OR M segment)",
        "(Rift Valley fever virus OR Phlebovirus riftense) AND (small segment OR S segment)"
    ]
    
    all_ids = set()  # Use set to avoid duplicates
    
    for strategy in search_strategies:
        print(f"\nSearching with strategy: {strategy}")
        try:
            # Get count first
            count_handle = Entrez.esearch(db="nucleotide", term=strategy, retmax=0)
            count_record = Entrez.read(count_handle)
            count_handle.close()
            
            total_count = int(count_record["Count"])
            print(f"Found {total_count} sequences for this search")
            
            if total_count == 0:
                continue
                
            # Get all IDs for this search
            actual_max = total_count if max_sequences is None else min(max_sequences, total_count)
            
            handle = Entrez.esearch(
                db="nucleotide", 
                term=strategy, 
                retmax=actual_max,
                usehistory="y"
            )
            record = Entrez.read(handle)
            handle.close()
            
            batch_ids = record["IdList"]
            print(f"Retrieved {len(batch_ids)} IDs")
            all_ids.update(batch_ids)
            
            time.sleep(0.5)  # Be nice to NCBI servers
            
        except Exception as e:
            print(f"Error with search strategy '{strategy}': {e}")
            continue
    
    final_ids = list(all_ids)
    print(f"\n=== SEARCH SUMMARY ===")
    print(f"Total unique sequences found: {len(final_ids)}")
    
    return final_ids

def fetch_sequences_with_full_metadata(id_list, batch_size=50):
    """
    Fetch sequences with comprehensive metadata matching NCBI Virus website format
    """
    sequences = []
    metadata_records = []
    segment_counts = {'L': 0, 'M': 0, 'S': 0, 'unknown': 0}
    
    print(f"\nFetching {len(id_list)} sequences with full metadata...")
    
    for start in range(0, len(id_list), batch_size):
        end = min(start + batch_size, len(id_list))
        batch_ids = id_list[start:end]
        
        print(f"Fetching batch {start//batch_size + 1}: {start+1}-{end} of {len(id_list)}")
        
        try:
            # Fetch GenBank records
            fetch_handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb", 
                retmode="text",
                id=",".join(batch_ids)
            )
            
            for record in SeqIO.parse(fetch_handle, "genbank"):
                # Initialize metadata with NCBI Virus website column structure
                metadata = {
                    "Accession": record.id,
                    "Organism_Name": "",
                    "GenBank_RefSeq": "GenBank",  # Default, will update based on accession
                    "Assembly": "",
                    "SRA_Accession": "",
                    "Submitters": "",
                    "Organization": "",
                    "Org_location": "",
                    "Release_Date": "",
                    "Isolate": "",
                    "Species": "",
                    "Genus": "",
                    "Family": "",
                    "Molecule_type": "",
                    "Length": len(record.seq),
                    "Nuc_Completeness": "",
                    "Genotype": "",
                    "Segment": "",
                    "Publications": "",
                    "Geo_Location": "",
                    "Country": "",
                    "USA": "",  # This column is typically empty
                    "Host": "",
                    "Tissue_Specimen_Source": "",
                    "Collection_Date": "",
                    "BioSample": "",
                    "BioProject": "",
                    "GenBank_Title": record.description
                }
                
                # Extract basic information
                if "date" in record.annotations:
                    metadata["Release_Date"] = record.annotations["date"]
                
                # Determine if RefSeq or GenBank
                if record.id.startswith(("NC_", "NG_", "NM_", "NP_", "NR_", "NT_", "NW_", "NZ_")):
                    metadata["GenBank_RefSeq"] = "RefSeq"
                
                # Extract information from dbxrefs
                if record.dbxrefs:
                    for xref in record.dbxrefs:
                        if xref.startswith("Assembly:"):
                            metadata["Assembly"] = xref.split(":", 1)[1]
                        elif xref.startswith("BioSample:"):
                            metadata["BioSample"] = xref.split(":", 1)[1]
                        elif xref.startswith("BioProject:"):
                            metadata["BioProject"] = xref.split(":", 1)[1]
                        elif xref.startswith("SRA:"):
                            metadata["SRA_Accession"] = xref.split(":", 1)[1]
                
                # Process features for detailed metadata
                for feature in record.features:
                    if feature.type == "source":
                        # Organism information
                        if "organism" in feature.qualifiers:
                            org_name = feature.qualifiers["organism"][0]
                            metadata["Organism_Name"] = org_name
                            
                            # Set correct species name based on current NCBI taxonomy
                            if "Rift Valley fever virus" in org_name or "Phlebovirus riftense" in org_name:
                                metadata["Species"] = "Phlebovirus riftense"  # Updated species name
                                metadata["Genus"] = "Phlebovirus"
                                metadata["Family"] = "Phenuiviridae"
                        
                        # Strain/Isolate
                        if "strain" in feature.qualifiers:
                            metadata["Isolate"] = feature.qualifiers["strain"][0]
                        elif "isolate" in feature.qualifiers:
                            metadata["Isolate"] = feature.qualifiers["isolate"][0]
                        
                        # Geographic information
                        if "country" in feature.qualifiers:
                            country_info = feature.qualifiers["country"][0]
                            metadata["Geo_Location"] = country_info
                            if ":" in country_info:
                                country, location = country_info.split(":", 1)
                                metadata["Country"] = country.strip()
                            else:
                                metadata["Country"] = country_info.strip()
                        
                        # Host information
                        if "host" in feature.qualifiers:
                            metadata["Host"] = feature.qualifiers["host"][0]
                        
                        # Collection date
                        if "collection_date" in feature.qualifiers:
                            metadata["Collection_Date"] = feature.qualifiers["collection_date"][0]
                        
                        # Tissue/specimen source
                        if "tissue_type" in feature.qualifiers:
                            metadata["Tissue_Specimen_Source"] = feature.qualifiers["tissue_type"][0]
                        elif "isolation_source" in feature.qualifiers:
                            metadata["Tissue_Specimen_Source"] = feature.qualifiers["isolation_source"][0]
                        
                        # Molecule type
                        if "mol_type" in feature.qualifiers:
                            metadata["Molecule_type"] = feature.qualifiers["mol_type"][0]
                
                # Determine completeness
                desc_lower = record.description.lower()
                if "complete genome" in desc_lower or "complete sequence" in desc_lower:
                    metadata["Nuc_Completeness"] = "complete"
                else:
                    metadata["Nuc_Completeness"] = "partial"
                
                # Segment identification - comprehensive approach
                segment_identified = False
                desc_lower = record.description.lower()
                
                # Check description for segment information
                if "segment l" in desc_lower or "large segment" in desc_lower:
                    metadata["Segment"] = "L"
                    segment_counts["L"] += 1
                    segment_identified = True
                elif "segment m" in desc_lower or "medium segment" in desc_lower or "middle segment" in desc_lower:
                    metadata["Segment"] = "M"
                    segment_counts["M"] += 1
                    segment_identified = True
                elif "segment s" in desc_lower or "small segment" in desc_lower:
                    metadata["Segment"] = "S"
                    segment_counts["S"] += 1
                    segment_identified = True
                
                # If not found in description, check features
                if not segment_identified:
                    for feature in record.features:
                        for qualifier_name in ["product", "note", "gene", "segment"]:
                            if qualifier_name in feature.qualifiers:
                                for qualifier_value in feature.qualifiers[qualifier_name]:
                                    qual_lower = qualifier_value.lower()
                                    if any(term in qual_lower for term in ["segment l", "large", "rdrp", "polymerase l"]):
                                        metadata["Segment"] = "L"
                                        segment_counts["L"] += 1
                                        segment_identified = True
                                        break
                                    elif any(term in qual_lower for term in ["segment m", "medium", "glyco", "gn", "gc"]):
                                        metadata["Segment"] = "M"
                                        segment_counts["M"] += 1
                                        segment_identified = True
                                        break
                                    elif any(term in qual_lower for term in ["segment s", "small", "nucleo", "n protein"]):
                                        metadata["Segment"] = "S"
                                        segment_counts["S"] += 1
                                        segment_identified = True
                                        break
                                if segment_identified:
                                    break
                        if segment_identified:
                            break
                
                # Default to unknown if still not identified
                if not segment_identified:
                    metadata["Segment"] = "unknown"
                    segment_counts["unknown"] += 1
                
                # Add to collections
                sequences.append(record)
                metadata_records.append(metadata)
            
            fetch_handle.close()
            time.sleep(0.5)  # Be nice to NCBI servers
            
        except Exception as e:
            print(f"Error fetching batch {start}-{end}: {e}")
            continue
    
    print(f"\n=== DOWNLOAD SUMMARY ===")
    print(f"Total sequences downloaded: {len(sequences)}")
    print(f"Segment distribution:")
    for segment, count in segment_counts.items():
        print(f"  {segment}: {count}")
    
    return sequences, metadata_records

def save_data(sequences, metadata_records, output_sequences, output_metadata):
    """Save sequences and metadata to files"""
    
    # Create output directories if they don't exist
    os.makedirs(os.path.dirname(output_sequences), exist_ok=True)
    os.makedirs(os.path.dirname(output_metadata), exist_ok=True)
    
    # Save sequences in FASTA format
    print(f"Writing {len(sequences)} sequences to {output_sequences}")
    with open(output_sequences, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    
    # Save metadata as TSV
    print(f"Writing metadata to {output_metadata}")
    df = pd.DataFrame(metadata_records)
    df.to_csv(output_metadata, sep="\t", index=False)
    
    return len(sequences), len(metadata_records)

def main():
    """Main function for standalone usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Download RVF sequences matching NCBI Virus website exactly")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--api-key", help="NCBI API key")
    parser.add_argument("--max-sequences", type=int, help="Maximum sequences to download")
    parser.add_argument("--output-sequences", default="data/sequences/raw/rvf_sequences.fasta")
    parser.add_argument("--output-metadata", default="data/metadata/raw/rvf_metadata.tsv")
    
    args = parser.parse_args()
    
    print("=== NCBI Virus Exact Data Download ===")
    print(f"Email: {args.email}")
    print(f"Max sequences: {args.max_sequences or 'All available'}")
    print(f"Output sequences: {args.output_sequences}")
    print(f"Output metadata: {args.output_metadata}")
    
    # Search for sequences
    id_list = search_ncbi_virus_exact(args.email, args.api_key, args.max_sequences)
    
    if not id_list:
        print("No sequences found!")
        return
    
    # Fetch sequences and metadata
    sequences, metadata_records = fetch_sequences_with_full_metadata(id_list)
    
    if not sequences:
        print("No sequences successfully downloaded!")
        return
    
    # Save data
    seq_count, meta_count = save_data(sequences, metadata_records, args.output_sequences, args.output_metadata)
    
    print(f"\n=== DOWNLOAD COMPLETE ===")
    print(f"Sequences saved: {seq_count}")
    print(f"Metadata records: {meta_count}")
    print(f"Files created:")
    print(f"  - {args.output_sequences}")
    print(f"  - {args.output_metadata}")

if __name__ == "__main__":
    main()
