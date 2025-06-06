#!/usr/bin/env python3
"""
Script to download RVF virus sequences and metadata from NCBI Virus
"""

import os
import argparse
import sys
import re
from datetime import datetime
from Bio import Entrez, SeqIO
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Download RVF sequences from NCBI Virus")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--api-key", help="NCBI API key for higher rate limits")
    parser.add_argument("--search-term", 
                        default="Rift Valley fever virus[organism] taxid:11588", 
                        help="Search term for NCBI Entrez")
    parser.add_argument("--max-sequences", type=int, default=None, 
                        help="Maximum number of sequences to download (default: all available)")
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

# Modified search_ncbi function with alternative search terms
def search_ncbi(search_term, max_results=None, email=None, api_key=None):
    """Search NCBI for sequences matching the search term"""
    if email:
        Entrez.email = email
    else:
        # This should not happen due to argparse 'required=True'
        print("Error: Email is required for NCBI Entrez.", file=sys.stderr)
        sys.exit(1) 
    
    if api_key:
        Entrez.api_key = api_key
        print(f"Using provided NCBI API key")
    search_terms = [
        search_term,  # Original search term first
        "Rift Valley fever virus[Organism] OR Phlebovirus riftense[Organism]",  # Include both organism names
        "txid11588[Organism:exp]",  # Direct taxid search
        "Rift Valley fever[Title] AND virus[Title]",  # Title-based search
        "Rift Valley fever virus OR Phlebovirus riftense",  # Plain text search including both names
        "Rift Valley fever"  # Even broader search
    ]
    
    for current_term in search_terms:
        print(f"Searching NCBI for: {current_term}")
        try:
            # First, get the total count of records without retrieving them
            count_handle = Entrez.esearch(db="nucleotide", term=current_term, retmax=0)
            count_record = Entrez.read(count_handle)
            count_handle.close()
            
            total_records = int(count_record["Count"])
            print(f"Found {total_records} total records matching the search criteria")
            
            if total_records == 0:
                if current_term == search_terms[-1]: # If this was the last term
                    print("All search terms failed to return results.")
                    return [], None, None
                print("Trying alternative search term...")
                continue  # Try next search term
                
            # Now retrieve all records (or up to max_results if specified)
            actual_max = total_records if max_results is None else min(max_results, total_records)
            
            handle = Entrez.esearch(db="nucleotide", term=current_term, 
                                  retmax=actual_max, usehistory="y")
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            web_env = record["WebEnv"]
            query_key = record["QueryKey"]
            
            print(f"Preparing to download {len(id_list)} records")
            return id_list, web_env, query_key
            
        except Exception as e:
            print(f"Error during NCBI search with term '{current_term}': {e}")
            if current_term == search_terms[-1]: # If this was the last term
                 print("All search terms failed.")
                 return [], None, None
            continue  # Try next search term
    
    # Should be unreachable if logic above is correct
    return [], None, None

def fetch_sequences(id_list, web_env=None, query_key=None, batch_size=100, debug=False):
    """Fetch sequences by ID list"""
    sequences = []
    metadata_records = []
    segment_counts = {'L': 0, 'M': 0, 'S': 0, 'unknown': 0}
    qualifier_values_debug_list = []
    
    if not id_list:
        print("No IDs to fetch.")
        return sequences, metadata_records, (qualifier_values_debug_list if debug else None)

    batch_size = min(batch_size, len(id_list)) if id_list else 1
    if batch_size == 0 and len(id_list) > 0:
        batch_size = 1

    # Define headers based on sequences.csv
    csv_headers = [
        "Accession", "Organism_Name", "GenBank_RefSeq", "Assembly", "SRA_Accession",
        "Submitters", "Organization", "Org_location", "Release_Date", "Isolate",
        "Species", "Genus", "Family", "Molecule_type", "Length", "Nuc_Completeness",
        "Genotype", "Segment", "Publications", "Geo_Location", "Country", "USA",
        "Host", "Tissue_Specimen_Source", "Collection_Date", "BioSample", "BioProject",
        "GenBank_Title"
    ]

    print(f"Fetching {len(id_list)} sequences")
    for start in range(0, len(id_list), batch_size):
        end = min(start + batch_size, len(id_list))
        batch_ids = id_list[start:end]
        
        print(f"Fetching batch {start//batch_size + 1} ({start+1}-{end}) of {len(id_list)}")
        
        try:
            if web_env and query_key:
                fetch_handle = Entrez.efetch(
                    db="nucleotide", rettype="gb", retmode="text",
                    retstart=start, retmax=batch_size, webenv=web_env, query_key=query_key
                )
            else:
                fetch_handle = Entrez.efetch(
                    db="nucleotide", rettype="gb", retmode="text", id=",".join(batch_ids)                )
                
            for record_idx, record in enumerate(SeqIO.parse(fetch_handle, "genbank")):
                metadata = {key: "" for key in csv_headers}

                metadata["Accession"] = record.id
                metadata["Length"] = len(record.seq)
                metadata["GenBank_Title"] = record.description
                metadata["USA"] = "" # USA column is always empty in the provided CSV

                if "date" in record.annotations:
                    metadata["Release_Date"] = record.annotations["date"]
                
                if "molecule_type" in record.annotations:
                    metadata["Molecule_type"] = record.annotations["molecule_type"]
                elif record.features:
                    for feature in record.features: # Check mol_type in features
                        if "mol_type" in feature.qualifiers:
                            metadata["Molecule_type"] = feature.qualifiers["mol_type"][0]
                            break
                
                desc_lower = record.description.lower()
                if "complete genome" in desc_lower:
                    metadata["Nuc_Completeness"] = "complete"
                else:
                    metadata["Nuc_Completeness"] = "partial"

                if record.id.startswith("NC_") or record.id.startswith("NG_") or \
                   record.id.startswith("NM_") or record.id.startswith("NP_") or \
                   record.id.startswith("NR_") or record.id.startswith("NT_") or \
                   record.id.startswith("NW_") or record.id.startswith("NZ_"):
                    metadata["GenBank_RefSeq"] = "RefSeq"
                else:
                    is_refseq_in_dbxrefs = False
                    if record.dbxrefs:
                        for xref in record.dbxrefs:
                            if "RefSeq" in xref: # A bit broad, but might catch some cases
                                is_refseq_in_dbxrefs = True
                                break
                    if is_refseq_in_dbxrefs:
                         metadata["GenBank_RefSeq"] = "RefSeq"
                    else:
                        metadata["GenBank_RefSeq"] = "GenBank"

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

                for feature in record.features:
                    if feature.type == "source":
                        if "organism" in feature.qualifiers:
                            org_name = feature.qualifiers["organism"][0]
                            metadata["Organism_Name"] = org_name
                            
                            # Standardize Species name regardless of how it appears in GenBank
                            if "Phlebovirus riftense" in org_name:
                                # For records using the newer name
                                metadata["Species"] = "Phlebovirus riftense"
                                # But add both names for compatibility
                                metadata["Organism_Name"] = "Rift Valley fever virus (Phlebovirus riftense)"
                            else:
                                # For records using the traditional name
                                metadata["Species"] = "Rift Valley fever virus"
                                
                            # Set genus and family consistently
                            metadata["Genus"] = "Phlebovirus"
                            metadata["Family"] = "Phenuiviridae"
                        
                        if "isolate" in feature.qualifiers:
                            metadata["Isolate"] = feature.qualifiers["isolate"][0]
                        elif "strain" in feature.qualifiers and not metadata["Isolate"]:
                            metadata["Isolate"] = feature.qualifiers["strain"][0]

                        if "collection_date" in feature.qualifiers:
                            metadata["Collection_Date"] = feature.qualifiers["collection_date"][0]
                        
                        if "country" in feature.qualifiers:
                            country_val = feature.qualifiers["country"][0]
                            metadata["Geo_Location"] = country_val
                            metadata["Country"] = country_val.split(":")[0].strip()
                        
                        if "host" in feature.qualifiers:
                            metadata["Host"] = feature.qualifiers["host"][0]
                        elif "lab_host" in feature.qualifiers:
                            metadata["Host"] = feature.qualifiers["lab_host"][0]
                            
                        if "tissue_type" in feature.qualifiers:
                            metadata["Tissue_Specimen_Source"] = feature.qualifiers["tissue_type"][0]
                        elif "isolation_source" in feature.qualifiers and not metadata["Tissue_Specimen_Source"]:
                            metadata["Tissue_Specimen_Source"] = feature.qualifiers["isolation_source"][0]

                        if "genotype" in feature.qualifiers:
                            metadata["Genotype"] = feature.qualifiers["genotype"][0]
                        elif "note" in feature.qualifiers:
                            for note_val in feature.qualifiers["note"]:
                                if "genotype" in note_val.lower():
                                    metadata["Genotype"] = note_val.split(":")[-1].strip() if ":" in note_val else note_val
                                    break
                        
                        if not metadata["Assembly"] and "db_xref" in feature.qualifiers:
                            for xref in feature.qualifiers["db_xref"]:
                                if xref.startswith("Assembly:"): metadata["Assembly"] = xref.split(":",1)[1]; break
                        if not metadata["BioSample"] and "db_xref" in feature.qualifiers:
                            for xref in feature.qualifiers["db_xref"]:
                                if xref.startswith("BioSample:"): metadata["BioSample"] = xref.split(":",1)[1]; break
                        if not metadata["BioProject"] and "db_xref" in feature.qualifiers:
                            for xref in feature.qualifiers["db_xref"]:
                                if xref.startswith("BioProject:"): metadata["BioProject"] = xref.split(":",1)[1]; break
                        if not metadata["SRA_Accession"] and "db_xref" in feature.qualifiers:
                            for xref in feature.qualifiers["db_xref"]:
                                if xref.startswith("SRA:"): metadata["SRA_Accession"] = xref.split(":",1)[1]; break
                        break 

                if record.annotations.get("references"):
                    first_ref = record.annotations["references"][0]
                    metadata["Submitters"] = first_ref.authors

                    journal_text = first_ref.journal
                    org_loc_parsed_from_journal = False
                    if "Submitted" in journal_text:
                        submission_info_match = re.search(r"Submitted\s*\([\d\w-]+\)\s*(.*)", journal_text, re.IGNORECASE)
                        if submission_info_match:
                            submission_details = submission_info_match.group(1).strip()
                            known_countries_map = { # Map for more robust matching, key is country name, value is part of address string
                                "USA": ["USA"], "Austria": ["Austria"], "Kenya": ["Kenya"], 
                                "Tanzania": ["Tanzania"], "Senegal": ["Senegal"], "Egypt": ["Egypt"],
                                "Botswana": ["Botswana"], "Rwanda": ["Rwanda"], "Uganda": ["Uganda"]
                            } # This can be expanded
                            
                            found_country_in_submission = ""
                            for country_name_map, country_strings_map in known_countries_map.items():
                                for country_str_map in country_strings_map:
                                    if submission_details.endswith(country_str_map):
                                        found_country_in_submission = country_name_map
                                        break
                                if found_country_in_submission:
                                    break
                            
                            if found_country_in_submission:
                                metadata["Org_location"] = found_country_in_submission
                                # Attempt to get text before the address part as organization
                                # This is a heuristic: find the last comma before a potential city name or long address string
                                # For simplicity, if country is found, take text before that as org.
                                # Example: "Org Name, Some Dept, City, Country" -> Org Name, Some Dept
                                # Example: "Org Name, Country" -> Org Name
                                # This needs to be robust. The CSV has very specific org names.
                                # Let's try to find the part of the string that is the org.
                                # The CSV for NC_014395.1 has "National Center for Biotechnology Information, NIH"
                                # The journal is "Submitted (11-AUG-2010) National Center for Biotechnology Information, NIH, Bethesda, MD 20894, USA"
                                # The CSV for PV231432.1 has "International Atomic Energy Agency, Animal Production and Health Laboratory, Joint FAO/IAEA Centre of Nuclear Techniques in Food and Agriculture, Department of Nuclear Sciences and Applications"
                                # The journal is "Submitted (23-MAR-2025) International Atomic Energy Agency, Animal Production and Health Laboratory, Joint FAO/IAEA Centre of Nuclear Techniques in Food and Agriculture, Department of Nuclear Sciences and Applications, Vienna International Centre, PO Box 100, Vienna, A-1400, Austria"
                                
                                # A simple split: find the first occurrence of a city-like pattern or the country itself if it's at the end of a segment.
                                # This is very hard to generalize perfectly.
                                # For now, a placeholder or a simpler split:
                                last_comma_before_country_guess = submission_details.rfind(',', 0, submission_details.lower().rfind(found_country_in_submission.lower()))
                                if last_comma_before_country_guess != -1:
                                     possible_org = submission_details[:last_comma_before_country_guess].strip()
                                     # Further refine: if org ends with a state/zip, try to remove that too.
                                     # e.g. "..., NIH, Bethesda, MD 20894" -> we want "..., NIH"
                                     state_zip_pattern = r",\s*[A-Za-z\s]+,\s*[A-Z]{2}\s*\d{5}(-\d{4})?$" # Matches ", City, ST ZIP"
                                     m_state_zip = re.search(state_zip_pattern, possible_org)
                                     if m_state_zip:
                                         possible_org = possible_org[:m_state_zip.start()]
                                     metadata["Organization"] = possible_org

                                else: # No comma before, maybe it's just "Org, Country"
                                    first_comma = submission_details.find(',')
                                    if first_comma != -1 and submission_details[first_comma+1:].strip().lower().startswith(found_country_in_submission.lower()):
                                        metadata["Organization"] = submission_details[:first_comma].strip()
                                    else:
                                         metadata["Organization"] = submission_details # Fallback
                                org_loc_parsed_from_journal = True

                            if not org_loc_parsed_from_journal: # Fallback if country parsing failed
                                metadata["Organization"] = submission_details
                    
                    has_publication = False
                    if record.annotations.get("references"):
                        for ref_item in record.annotations["references"]:
                            if (hasattr(ref_item, 'pubmed_id') and ref_item.pubmed_id and ref_item.pubmed_id != "0") or \
                               (hasattr(ref_item, 'title') and "direct submission" not in ref_item.title.lower()):
                                has_publication = True
                                break
                    if has_publication:
                        metadata["Publications"] = "1"
                  # Segment determination (enhanced detection for RVF virus)
                current_segment = ""
                segment_found = False
                
                # 1. Look for segment in GenBank "segment" feature qualifier
                for feature in record.features:
                    if "segment" in feature.qualifiers:
                        segment_value = feature.qualifiers["segment"][0].upper()
                        if segment_value in ["L", "M", "S"]:
                            current_segment = segment_value
                            segment_found = True
                            break
                
                # 2. From CDS product/gene names if segment qualifier not found
                if not segment_found:
                    for feature in record.features:
                        for qualifier_key in ["product", "gene", "note"]:
                            if qualifier_key in feature.qualifiers:
                                for value in feature.qualifiers[qualifier_key]:
                                    value_lower = value.lower()
                                    
                                    # L segment detection patterns
                                    if any(k in value_lower for k in ["rdrp", "rna-dependent rna polymerase", "l protein", 
                                                                    "polymerase", "large", "l segment", "segment l"]):
                                        current_segment = "L"
                                        segment_found = True
                                        break
                                    
                                    # M segment detection patterns
                                    elif any(k in value_lower for k in ["glycoprotein", "gc", "gn", "m protein", 
                                                                      "envelope protein", "medium", "m segment", "segment m"]):
                                        current_segment = "M"
                                        segment_found = True
                                        break
                                    
                                    # S segment detection patterns
                                    elif any(k in value_lower for k in ["nucleocapsid", "np", "n protein", "nss", 
                                                                      "non-structural protein s", "small", "s segment", "segment s"]):
                                        current_segment = "S"
                                        segment_found = True
                                        break
                                if segment_found: break
                        if segment_found: break
                
                # 3. From record description if not found yet
                if not segment_found and record.description:
                    description_lower = record.description.lower()
                    
                    # Check for segment identifiers in the description
                    if any(pattern in description_lower for pattern in ["segment l", "l segment", "l rna", "large segment", "l protein gene"]):
                        current_segment = "L"
                        segment_found = True
                    elif any(pattern in description_lower for pattern in ["segment m", "m segment", "m rna", "medium segment", "glycoprotein gene"]):
                        current_segment = "M"
                        segment_found = True
                    elif any(pattern in description_lower for pattern in ["segment s", "s segment", "s rna", "small segment", "nucleocapsid gene"]):
                        current_segment = "S"
                        segment_found = True
                
                # 4. Fallback to length-based identification if still not found
                if not segment_found:
                    seq_len = len(record.seq)
                    if seq_len > 5500:  # L segment typical length
                        current_segment = "L"
                    elif seq_len > 3000:  # M segment typical length
                        current_segment = "M"
                    elif seq_len > 800:  # S segment typical length
                        current_segment = "S"
                    else:
                        current_segment = "unknown"
                
                metadata["Segment"] = current_segment
                segment_counts[current_segment] = segment_counts.get(current_segment, 0) + 1

                sequences.append(record)
                metadata_records.append(metadata)
                
                if debug:
                    all_qualifier_values_for_record = []
                    for feature_debug in record.features:
                        for qual_key, qual_vals in feature_debug.qualifiers.items():
                            all_qualifier_values_for_record.extend([f"{qual_key}: {qv}" for qv in qual_vals])
                    qualifier_values_debug_list.append({
                        "accession": metadata["Accession"],
                        "segment_identified": metadata["Segment"],
                        "length": metadata["Length"],
                        "description": metadata["GenBank_Title"],
                        "completeness_identified": metadata["Nuc_Completeness"],
                        "qualifiers_sample": all_qualifier_values_for_record[:10]
                    })
                
            fetch_handle.close()
            
        except Exception as e:
            print(f"Error fetching or parsing batch {start+1}-{end}: {e}")
            # import traceback
            # traceback.print_exc() # For detailed error during debugging
            continue 
    
    print(f"Segment counts: L={segment_counts['L']}, M={segment_counts['M']}, "
          f"S={segment_counts['S']}, unknown={segment_counts['unknown']}")
    
    if debug and qualifier_values_debug_list:
        print("\nDebug: Sample of qualifier values from records (first 5):")
        for i, record_data_debug in enumerate(qualifier_values_debug_list[:5]): # Renamed record_data
            print(f"\nAccession: {record_data_debug['accession']}")
            print(f"Identified segment: {record_data_debug['segment_identified']}")
            print(f"Identified completeness: {record_data_debug['completeness_identified']}")
            print(f"Sequence length: {record_data_debug['length']}")
            print(f"Description: {record_data_debug['description']}")
            print("Sample qualifier values:")
            # Check if 'qualifiers_sample' exists and is not empty
            if record_data_debug.get('qualifiers_sample'):
                for val in record_data_debug['qualifiers_sample']:
                    print(f"  {val}")
                # Check if all_qualifier_values_for_record is defined in this scope for comparison
                # This comparison might be problematic if all_qualifier_values_for_record is not available here
                # For simplicity, just indicate if there were more than sampled.
                if len(record_data_debug['qualifiers_sample']) == 10 : print("  ...") 
            else: 
                print("  (No qualifiers found or sampled for this record)")
    
    return sequences, metadata_records, (qualifier_values_debug_list if debug else None)

def create_dummy_data(segment_choice, output_sequences, output_metadata, headers):
    """Create dummy data for testing when no real sequences are found"""
    segment = segment_choice if segment_choice != 'all' else 'L'
    print(f"Creating dummy {segment} segment data for pipeline testing")
    
    seq_length, nuc_completeness_val = 1000, "partial"
    if segment == 'L': seq_length, nuc_completeness_val = 6404, "complete" # Assuming L is whole genome
    elif segment == 'M': seq_length, nuc_completeness_val = 3885, "partial"
    elif segment == 'S': seq_length, nuc_completeness_val = 1690, "partial"
    
    with open(output_sequences, 'w') as f:
        f.write(f">DUMMY_{segment}_ACC|Rift Valley fever virus dummy {segment} segment|{nuc_completeness_val}\n")
        f.write("A" * seq_length + "\n")
    
    with open(output_metadata, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        dummy_row = {key: "" for key in headers}
        dummy_row.update({
            "Accession": f"DUMMY_{segment}_ACC",
            "Organism_Name": "Rift Valley fever virus",
            "GenBank_RefSeq": "GenBank",
            "Assembly": "DUMMY_ASM",
            "Submitters": "Dummy Submitter",
            "Organization": "Dummy Organization",
            "Org_location": "DummyLocation",
            "Release_Date": datetime.now().strftime("%Y-%m-%d"),
            "Isolate": f"RVFV_dummy_{segment}_isolate",
            "Species": "Phlebovirus riftense",
            "Genus": "Phlebovirus",
            "Family": "Phenuiviridae",
            "Molecule_type": "ssRNA",
            "Length": str(seq_length),
            "Nuc_Completeness": nuc_completeness_val,
            "Segment": segment,
            "Geo_Location": "Kenya",
            "Country": "Kenya",
            "Host": "Bos taurus",
            "Collection_Date": datetime.now().strftime("%Y-%m-%d"),
            "BioSample": "DUMMY_BIOSAMPLE",
            "BioProject": "DUMMY_BIOPROJECT",
            "GenBank_Title": f"Rift Valley fever virus dummy {segment} segment, {nuc_completeness_val} sequence"
        })
        f.write('\t'.join([str(dummy_row.get(h, "")) for h in headers]) + '\n')

def main():
    args = parse_args()
    
    os.makedirs(os.path.dirname(args.output_sequences), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_metadata), exist_ok=True)
        
    id_list, web_env, query_key = search_ncbi(
        args.search_term, args.max_sequences, args.email, args.api_key
    )
    
    # Headers based on sequences.csv
    csv_headers = [
        "Accession", "Organism_Name", "GenBank_RefSeq", "Assembly", "SRA_Accession",
        "Submitters", "Organization", "Org_location", "Release_Date", "Isolate",
        "Species", "Genus", "Family", "Molecule_type", "Length", "Nuc_Completeness",
        "Genotype", "Segment", "Publications", "Geo_Location", "Country", "USA",
        "Host", "Tissue_Specimen_Source", "Collection_Date", "BioSample", "BioProject",
        "GenBank_Title"
    ]

    if not id_list:
        print("No sequences found matching the search term.")
        if args.create_dummy:
            create_dummy_data(args.segment, args.output_sequences, args.output_metadata, csv_headers)
        else: 
            with open(args.output_sequences, 'w') as f: pass
            with open(args.output_metadata, 'w') as f:
                f.write('\t'.join(csv_headers) + '\n')
        return 0 
    
    sequences, metadata_records, _ = fetch_sequences(
        id_list, web_env, query_key, debug=args.debug
    )
    
    if not sequences:
        print("Failed to fetch any sequences from the identified IDs.")
        if args.create_dummy:
            create_dummy_data(args.segment, args.output_sequences, args.output_metadata, csv_headers)
        else: 
            with open(args.output_sequences, 'w') as f: pass
            with open(args.output_metadata, 'w') as f:
                f.write('\t'.join(csv_headers) + '\n')
        return 0 
    
    if args.segment != 'all':
        filtered_sequences = [
            seq for seq, meta in zip(sequences, metadata_records) 
            if meta.get("Segment") == args.segment # Updated key
        ]
        filtered_metadata = [
            meta for meta in metadata_records 
            if meta.get("Segment") == args.segment # Updated key
        ]
        
        sequences, metadata_records = filtered_sequences, filtered_metadata
        print(f"Filtered to {len(sequences)} sequences of segment {args.segment}")
    
    if not sequences: 
        print(f"No sequences found for segment {args.segment} after filtering (if applicable).")
        if args.create_dummy:
            create_dummy_data(args.segment, args.output_sequences, args.output_metadata, csv_headers)
        else: 
            with open(args.output_sequences, 'w') as f: pass
            with open(args.output_metadata, 'w') as f:
                f.write('\t'.join(csv_headers) + '\n')
        return 0

    with open(args.output_sequences, 'w') as f:
        SeqIO.write(sequences, f, 'fasta')
    
    if metadata_records:
        df = pd.DataFrame(metadata_records)
        df = df.reindex(columns=csv_headers) 
        df.to_csv(args.output_metadata, sep='\t', index=False, na_rep='') # Use empty string for NA as in CSV
    else: 
        with open(args.output_metadata, 'w') as f:
            f.write('\t'.join(csv_headers) + '\n')
    
    print(f"Downloaded {len(sequences)} sequences to {args.output_sequences}")
    print(f"Wrote metadata for {len(metadata_records)} sequences to {args.output_metadata}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())