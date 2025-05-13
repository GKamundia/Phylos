#!/usr/bin/env python3
"""
Script to download RVF virus sequences and metadata from NCBI
"""
import os
import sys
import argparse
import time
import urllib.error
import http.client
import socket
from datetime import datetime
from Bio import Entrez, SeqIO
import pandas as pd

# Configure basic logging
import logging
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def retry(func, max_tries=5, delay=2.0, backoff=2.0, exceptions=(Exception,)):
    """Simple retry decorator for API calls"""
    def wrapper(*args, **kwargs):
        mtries, mdelay = max_tries, delay
        last_exception = None
        
        while mtries > 0:
            try:
                return func(*args, **kwargs)
            except exceptions as e:
                last_exception = e
                logger.warning(f"Request failed: {str(e)}, retrying in {mdelay}s... ({mtries} tries left)")
                time.sleep(mdelay)
                mtries -= 1
                mdelay *= backoff
        
        raise last_exception
    
    return wrapper

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

@retry
def search_ncbi(search_term, max_results=1000, email=None, api_key=None):
    """Search NCBI for sequences matching the search term with retry logic"""
    # Set Entrez email and API key
    if email:
        Entrez.email = email
    else:
        Entrez.email = "example@example.com"
    
    # Set API key if provided
    if api_key:
        Entrez.api_key = api_key
        logger.info("Using provided NCBI API key")
    
    logger.info(f"Searching NCBI for: {search_term}")
    
    try:
        # First just attempt to get the count to see if search backend is working
        count_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=5)
        count_record = Entrez.read(count_handle)
        count_handle.close()
        
        total_count = int(count_record["Count"])
        logger.info(f"Found {total_count} sequences total")
        
        # Now get the actual results with WebEnv
        handle = Entrez.esearch(
            db="nucleotide", 
            term=search_term, 
            retmax=max_results, 
            usehistory="y"
        )
        
        try:
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            web_env = record["WebEnv"]
            query_key = record["QueryKey"]
            
            logger.info(f"Retrieved {len(id_list)} sequence IDs")
            return id_list, web_env, query_key
            
        except RuntimeError as e:
            logger.error(f"NCBI search backend error: {str(e)}")
            # If we got the count earlier, we can try a different approach
            if total_count > 0:
                logger.info("Trying to fetch by batch without WebEnv...")
                # Try direct ID retrieval without WebEnv
                batch_size = 200
                id_list = []
                for start in range(0, min(total_count, max_results), batch_size):
                    logger.info(f"Fetching batch starting at {start}")
                    batch_handle = Entrez.esearch(
                        db="nucleotide", 
                        term=search_term, 
                        retstart=start, 
                        retmax=min(batch_size, max_results - len(id_list))
                    )
                    batch_record = Entrez.read(batch_handle)
                    batch_handle.close()
                    id_list.extend(batch_record["IdList"])
                    if len(id_list) >= max_results:
                        break
                    time.sleep(1)  # Be gentle to the NCBI servers
                
                return id_list, None, None
            else:
                raise
            
    except (urllib.error.URLError, http.client.HTTPException, socket.timeout) as e:
        logger.error(f"Network error: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during NCBI search: {str(e)}")
        raise

def fetch_sequences(id_list, web_env=None, query_key=None, batch_size=100):
    """Fetch sequences by ID list"""
    sequences = []
    metadata_records = []
    segment_counts = {'L': 0, 'M': 0, 'S': 0, 'unknown': 0}
    batch_size = min(batch_size, len(id_list))
    
    logger.info(f"Fetching {len(id_list)} sequences")
    for start in range(0, len(id_list), batch_size):
        end = min(start + batch_size, len(id_list))
        batch_ids = id_list[start:end]
        
        logger.info(f"Fetching batch {start//batch_size + 1} ({start}-{end})")
        
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
            
            # Sleep between batches to avoid overwhelming NCBI
            if end < len(id_list):
                time.sleep(1)
                
        except Exception as e:
            logger.error(f"Error fetching batch {start}-{end}: {e}")
            continue
    
    logger.info(f"Segment counts: L={segment_counts['L']}, M={segment_counts['M']}, "
          f"S={segment_counts['S']}, unknown={segment_counts['unknown']}")
    
    return sequences, metadata_records

def main():
    args = parse_args()
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(args.output_sequences), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_metadata), exist_ok=True)
    
    try:
        # Search NCBI with retry logic
        logger.info("Starting NCBI search...")
        id_list, web_env, query_key = search_ncbi(
            args.search_term, 
            args.max_sequences, 
            args.email,
            args.api_key
        )
        
        if not id_list:
            logger.error("No sequences found matching the search term")
            # Create empty output files so workflow can continue
            with open(args.output_sequences, 'w') as f:
                f.write(">dummy_sequence\nACGT\n")
            
            with open(args.output_metadata, 'w') as f:
                f.write("strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength\n")
                f.write("dummy_sequence\tRift Valley fever virus\tDUMMY1\t2023-01-01\tUnknown\t\t\t\tL\t4\n")
            
            logger.info(f"Created empty placeholder files to allow pipeline to continue")
            return 1
        
        # Fetch sequences and metadata
        logger.info("Fetching sequence data...")
        sequences, metadata_records = fetch_sequences(id_list, web_env, query_key)
        
        if not sequences:
            logger.error("Failed to fetch any sequences")
            # Create empty output files so workflow can continue
            with open(args.output_sequences, 'w') as f:
                f.write(">dummy_sequence\nACGT\n")
            
            with open(args.output_metadata, 'w') as f:
                f.write("strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength\n")
                f.write("dummy_sequence\tRift Valley fever virus\tDUMMY1\t2023-01-01\tUnknown\t\t\t\tL\t4\n")
            
            logger.info(f"Created empty placeholder files to allow pipeline to continue")
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
            logger.info(f"Filtered to {len(sequences)} sequences of segment {args.segment}")
        
        # If still no sequences after filtering, create placeholder
        if not sequences:
            logger.warning(f"No sequences found for segment {args.segment}")
            # Create empty output files so workflow can continue
            with open(args.output_sequences, 'w') as f:
                f.write(">dummy_sequence\nACGT\n")
            
            with open(args.output_metadata, 'w') as f:
                f.write("strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength\n")
                f.write(f"dummy_sequence\tRift Valley fever virus\tDUMMY1\t2023-01-01\tUnknown\t\t\t\t{args.segment}\t4\n")
            
            logger.info(f"Created segment placeholder files to allow pipeline to continue")
            return 1
        
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
        
        logger.info(f"Downloaded {len(sequences)} sequences to {args.output_sequences}")
        logger.info(f"Wrote metadata to {args.output_metadata}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Unhandled exception: {e}")
        
        # Create fallback files to allow pipeline to continue
        logger.info("Creating fallback files to allow pipeline to continue")
        with open(args.output_sequences, 'w') as f:
            f.write(">dummy_sequence\nACGT\n")
        
        with open(args.output_metadata, 'w') as f:
            f.write("strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength\n")
            f.write(f"dummy_sequence\tRift Valley fever virus\tDUMMY1\t2023-01-01\tUnknown\t\t\t\t{args.segment if args.segment != 'all' else 'L'}\t4\n")
        
        return 1

if __name__ == "__main__":
    sys.exit(main())