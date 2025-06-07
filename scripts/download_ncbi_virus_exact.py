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
    else:
        # Use the new provided API key by default
        Entrez.api_key = "f3edafc1dc58cc4e438d90134e9191616b08"
        print(f"Using new NCBI API key for enhanced access")

def search_ncbi_virus_exact(email, api_key=None, max_sequences=None):
    """
    Search NCBI using the exact same approach as the NCBI Virus website
    This ensures we get ALL RVF virus sequences including all segments
    """
    setup_entrez(email, api_key)
    
    # Use the exact same search strategy as NCBI Virus website
    # This matches the website query: txid11588[Organism:exp] which gets ALL RVF sequences
    search_strategies = [
        # Primary search - EXACT match to NCBI Virus website
        "txid11588[Organism:exp]",  # This is the exact query used by NCBI Virus website
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
    Fetch sequences with comprehensive metadata directly from NCBI without any manipulation
    """
    sequences = []
    metadata_records = []
    
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
                            
                            # Set species/genus/family information directly from taxonomy
                            if "Rift Valley fever virus" in org_name or "Phlebovirus riftense" in org_name:
                                metadata["Species"] = "Phlebovirus riftense"
                                metadata["Genus"] = "Phlebovirus"
                                metadata["Family"] = "Phenuiviridae"
                        
                        # Strain/Isolate
                        if "strain" in feature.qualifiers:
                            metadata["Isolate"] = feature.qualifiers["strain"][0]
                        elif "isolate" in feature.qualifiers:
                            metadata["Isolate"] = feature.qualifiers["isolate"][0]
                        
                        # Improved Geographic information
                        country_info = feature.qualifiers.get("country", [None])[0]
                        if not country_info: # Fallback
                            country_info = feature.qualifiers.get("geo_loc_name", [None])[0]
                        
                        if country_info:
                            metadata["Geo_Location"] = country_info # Store the full string
                            if ":" in country_info:
                                country_part, 나머지는_무시 = country_info.split(":", 1)
                                metadata["Country"] = country_part.strip()
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
                        
                        # Molecule type - use raw value directly from NCBI
                        if "mol_type" in feature.qualifiers:
                            metadata["Molecule_type"] = feature.qualifiers["mol_type"][0]

                        # Segment from feature qualifier (preferred)
                        segment_qual = feature.qualifiers.get("segment", [None])[0]
                        if segment_qual:
                            # Ensure consistent casing (e.g., S, M, L)
                            metadata["Segment"] = segment_qual.upper().replace("SEGMENT ", "")

                # Extract Submitter and Organization from References
                if record.annotations.get("references"):
                    primary_ref = None
                    direct_submission_ref = None
                    
                    for ref_obj in record.annotations["references"]:
                        if not primary_ref: 
                            primary_ref = ref_obj
                        if hasattr(ref_obj, 'title') and "Direct Submission" in ref_obj.title:
                            direct_submission_ref = ref_obj
                            break 
                    
                    chosen_ref = direct_submission_ref if direct_submission_ref else primary_ref

                    if chosen_ref:
                        if hasattr(chosen_ref, 'authors') and chosen_ref.authors:
                             metadata["Submitters"] = chosen_ref.authors
                        
                        if hasattr(chosen_ref, 'journal') and chosen_ref.journal:
                            # For "Direct Submission", affiliation is often in the journal line:
                            # JOURNAL Submitted (DD-MON-YYYY) Department of Virology, Institute of Research...
                            journal_text = chosen_ref.journal
                            if "Submitted" in journal_text:
                                match = re.search(r'Submitted\s*\(.+?\)\s*(.*)', journal_text)
                                if match and match.group(1):
                                    metadata["Organization"] = match.group(1).strip()
                        
                        # Fallback or alternative for Organization if not found in journal
                        if not metadata["Organization"] and hasattr(chosen_ref, 'consrtm') and chosen_ref.consrtm:
                            metadata["Organization"] = chosen_ref.consrtm
                        elif not metadata["Organization"] and hasattr(chosen_ref, 'comment') and chosen_ref.comment:
                            # The comment field can sometimes hold affiliation details.
                            # This is less structured, so use with caution or add more parsing.
                            # For now, if Organization is still empty, we can try this.
                            # Example: "Protein Science Center,Tsinghua University"
                            # Avoid overly long or generic comments.
                            if len(chosen_ref.comment) < 200 and ("university" in chosen_ref.comment.lower() or "institute" in chosen_ref.comment.lower()):
                                metadata["Organization"] = chosen_ref.comment.strip()

                # Extract completeness information from GenBank record
                # This should match exactly what NCBI Virus website shows
                title_lower = record.description.lower()
                if 'complete genome' in title_lower or 'complete sequence' in title_lower:
                    metadata["Nuc_Completeness"] = "complete"
                elif 'partial' in title_lower:
                    metadata["Nuc_Completeness"] = "partial"
                else:
                    # Check annotations for completeness
                    if "completeness" in record.annotations:
                        comp = record.annotations["completeness"].lower()
                        if "complete" in comp:
                            metadata["Nuc_Completeness"] = "complete"
                        else:
                            metadata["Nuc_Completeness"] = "partial"
                    else:
                        # Default to partial if unclear
                        metadata["Nuc_Completeness"] = "partial"
                
                # Extract segment information from description (fallback if not in feature)
                if not metadata["Segment"]: # If not found in feature qualifier
                    title_lower = record.description.lower()
                    if 'segment s' in title_lower or 's segment' in title_lower:
                        metadata["Segment"] = "S"
                    elif 'segment m' in title_lower or 'm segment' in title_lower:
                        metadata["Segment"] = "M"
                    elif 'segment l' in title_lower or 'l segment' in title_lower:
                        metadata["Segment"] = "L"
                
                # NO DATA MANIPULATION - just store raw sequence and metadata
                sequences.append(record)
                metadata_records.append(metadata)
            
            fetch_handle.close()
            time.sleep(0.3)  # Be nice to NCBI
            
        except Exception as e:
            print(f"Error fetching batch {start//batch_size + 1}: {e}")
            continue
    
    return sequences, metadata_records

def save_data(sequences, metadata_records, output_fasta, output_metadata):
    """Save sequences and metadata to files"""
    print(f"\nSaving {len(sequences)} sequences to {output_fasta}")
    with open(output_fasta, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    
    print(f"Saving {len(metadata_records)} metadata records to {output_metadata}")
    df = pd.DataFrame(metadata_records)
    df.to_csv(output_metadata, sep='\t', index=False)
    
    print("\n=== DOWNLOAD COMPLETE ===")
    print(f"Total sequences downloaded: {len(sequences)}")
    print(f"FASTA file: {output_fasta}")
    print(f"Metadata file: {output_metadata}")

def main():
    """Main function to run the download"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Download RVF virus sequences from NCBI")
    parser.add_argument("--email", required=True, help="Your email address for NCBI")
    parser.add_argument("--api-key", help="NCBI API key (optional but recommended)")
    parser.add_argument("--max-sequences", type=int, help="Maximum number of sequences to download")
    parser.add_argument("--output-fasta", default="data/sequences/raw/rvf_sequences.fasta", 
                       help="Output FASTA file path")
    parser.add_argument("--output-metadata", default="data/metadata/raw/rvf_metadata.tsv",
                       help="Output metadata file path")
    
    args = parser.parse_args()
    
    # Create output directories if they don't exist
    os.makedirs(os.path.dirname(args.output_fasta), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_metadata), exist_ok=True)
    
    print("=== RVF VIRUS SEQUENCE DOWNLOAD ===")
    print(f"Using exact NCBI Virus website search: txid11588[Organism:exp]")
    print(f"Email: {args.email}")
    if args.api_key:
        print("Using API key for enhanced access")
    
    # Search for sequences
    sequence_ids = search_ncbi_virus_exact(
        email=args.email,
        api_key=args.api_key,
        max_sequences=args.max_sequences
    )
    
    if not sequence_ids:
        print("No sequences found!")
        return
    
    # Fetch sequences and metadata
    sequences, metadata_records = fetch_sequences_with_full_metadata(sequence_ids)
    
    # Save data
    save_data(sequences, metadata_records, args.output_fasta, args.output_metadata)

if __name__ == "__main__":
    main()
