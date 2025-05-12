#!/usr/bin/env python3
"""
Enhanced sequence download utility for pathogen genomic surveillance
Supports configurable data sources, versioned archiving, and incremental updates
"""

import logging
import os
import sys
import json
import time
import argparse
import datetime
from pathlib import Path
import shutil
from typing import Dict, List, Optional, Union
from Bio import Entrez
from Bio import SeqIO

# Import our new logging utilities
sys.path.append(os.path.join(os.path.dirname(__file__), '.'))
from utils.log_utils import setup_logger, log_execution_stats, log_with_context

# Configure logger
logger = setup_logger(
    name="download_sequences",
    log_file="download_sequences.log",
    level="INFO",
    json_output=True
)

class DataTracker:
    """Tracks downloaded data for incremental updates"""
    
    def __init__(self, tracking_file):
        self.tracking_file = tracking_file
        self.data = self._load_tracking_data()
    
    def _load_tracking_data(self) -> Dict:
        """Load tracking data from file or initialize if it doesn't exist"""
        try:
            if os.path.exists(self.tracking_file):
                with open(self.tracking_file, "r") as f:
                    data = json.load(f)
                log_with_context(logger, "INFO", f"Loaded tracking data from {self.tracking_file}", 
                                {"records": len(data.get("sequences", {})), "last_update": data.get("last_update")})
                return data
            else:
                log_with_context(logger, "INFO", f"No tracking file found at {self.tracking_file}, initializing new data")
                return {
                    "last_update": None,
                    "sequences": {"total_count": 0, "last_sequence_id": None},
                    "segments": {
                        "L": {"count": 0, "last_update": None},
                        "M": {"count": 0, "last_update": None},
                        "S": {"count": 0, "last_update": None}
                    }
                }
        except Exception as e:
            log_with_context(logger, "ERROR", f"Error loading tracking data: {e}", {"file": self.tracking_file})
            return {
                "last_update": None,
                "sequences": {"total_count": 0, "last_sequence_id": None},
                "segments": {
                    "L": {"count": 0, "last_update": None},
                    "M": {"count": 0, "last_update": None},
                    "S": {"count": 0, "last_update": None}
                }
            }
    
    def save(self) -> None:
        """Save tracking data to file"""
        try:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.tracking_file), exist_ok=True)
            
            with open(self.tracking_file, "w") as f:
                json.dump(self.data, f, indent=2)
            
            log_with_context(logger, "INFO", f"Saved tracking data to {self.tracking_file}", 
                           {"sequences": self.data["sequences"]["total_count"]})
        except Exception as e:
            log_with_context(logger, "ERROR", f"Failed to save tracking data: {e}", {"file": self.tracking_file})
    
    def update(self, stats: Dict) -> None:
        """Update tracking data with new download stats"""
        # Update last update timestamp
        self.data["last_update"] = datetime.datetime.now().strftime("%Y-%m-%d")
        
        # Update sequence counts
        self.data["sequences"]["total_count"] += stats.get("total", 0)
        
        # Update last sequence ID if provided
        if stats.get("last_id"):
            self.data["sequences"]["last_sequence_id"] = stats["last_id"]
        
        # Update segment counts if provided
        for segment, count in stats.get("segments", {}).items():
            if segment in self.data["segments"]:
                self.data["segments"][segment]["count"] += count
                self.data["segments"][segment]["last_update"] = self.data["last_update"]
        
        self.save()
    
    def get_last_update_date(self) -> Optional[str]:
        """Get the date of the last update"""
        return self.data.get("last_update")

class NCBIDataSource:
    """Handles sequence fetching from NCBI GenBank"""
    
    def __init__(self, email: str):
        self.email = email
        Entrez.email = email
        log_with_context(logger, "INFO", f"Initialized NCBI data source with email: {email}")
    
    def _search_sequences(self, query: str, max_results: int = 500, from_date: Optional[str] = None) -> List[str]:
        """
        Search for sequence IDs matching the query
        
        Args:
            query: Search term for NCBI
            max_results: Maximum number of results to return
            from_date: Only return results from this date forward (YYYY-MM-DD)
            
        Returns:
            List of sequence IDs
        """
        try:
            # Add date constraint if provided
            if from_date:
                date_obj = datetime.datetime.strptime(from_date, "%Y-%m-%d")
                # Format date for NCBI query: YYYY/MM/DD
                ncbi_date = date_obj.strftime("%Y/%m/%d")
                query = f"{query} AND {ncbi_date}[PDAT] : 3000[PDAT]"
            
            log_with_context(logger, "INFO", f"Searching NCBI with query: {query}", 
                           {"max_results": max_results, "from_date": from_date})
            
            # Search for records
            search_handle = Entrez.esearch(
                db="nucleotide",
                term=query,
                retmax=max_results,
                usehistory="y"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            # Extract IDs
            id_list = search_results["IdList"]
            
            # Log search stats
            log_with_context(logger, "INFO", f"Found {len(id_list)} sequences matching query", 
                           {"total_found": search_results["Count"], "returned": len(id_list)})
            
            return id_list, search_results.get("WebEnv"), search_results.get("QueryKey")
        
        except Exception as e:
            log_with_context(logger, "ERROR", f"Error searching NCBI: {e}", {"query": query})
            return [], None, None
    
    def _fetch_sequences(self, id_list: List[str], web_env=None, query_key=None, segment=None) -> Dict:
        """Fetch sequences by ID list or history"""
        if not id_list and not (web_env and query_key):
            log_with_context(logger, "WARNING", "No sequence IDs to fetch and no history provided")
            return {"sequences": [], "metadata_lines": ["strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength"], "stats": {"total": 0, "segments": {}}}
        
        try:
            batch_size = 100
            sequences = []
            metadata_lines = ["strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength"]
            segment_counts = {"L": 0, "M": 0, "S": 0, "unknown": 0}
            last_id = None
            
            log_with_context(logger, "INFO", f"Fetching {len(id_list)} sequences from NCBI", 
                           {"batch_size": batch_size})
            
            # Process in batches to avoid overloading NCBI's servers
            for start in range(0, len(id_list), batch_size):
                end = min(start + batch_size, len(id_list))
                batch_ids = id_list[start:end]
                
                # Log batch progress
                log_with_context(logger, "INFO", f"Fetching batch {start//batch_size + 1} ({len(batch_ids)} sequences)")
                
                # Fetch this batch
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
                
                # Parse records
                records = list(SeqIO.parse(fetch_handle, "genbank"))
                fetch_handle.close()
                
                # Process each record
                for record in records:
                    # Extract metadata
                    accession = record.id
                    last_id = accession
                    
                    # Basic metadata fields with defaults
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
                        "length": len(record)
                    }
                    
                    # Extract collection date
                    for feature in record.features:
                        if feature.type == "source":
                            # Date
                            if "collection_date" in feature.qualifiers:
                                metadata["date"] = feature.qualifiers["collection_date"][0]
                            
                            # Country
                            if "country" in feature.qualifiers:
                                metadata["country"] = feature.qualifiers["country"][0]
                            
                            # Host
                            if "host" in feature.qualifiers:
                                metadata["host"] = feature.qualifiers["host"][0]
                            
                            # Try to determine segment from product or note qualifiers
                            segment_found = False
                            for qualifier in ["product", "note", "segment"]:
                                if qualifier in feature.qualifiers:
                                    qualifier_value = feature.qualifiers[qualifier][0].lower()
                                    if "segment l" in qualifier_value or "large" in qualifier_value:
                                        metadata["segment"] = "L"
                                        segment_counts["L"] += 1
                                        segment_found = True
                                        break
                                    elif "segment m" in qualifier_value or "medium" in qualifier_value:
                                        metadata["segment"] = "M"
                                        segment_counts["M"] += 1
                                        segment_found = True
                                        break
                                    elif "segment s" in qualifier_value or "small" in qualifier_value:
                                        metadata["segment"] = "S"
                                        segment_counts["S"] += 1
                                        segment_found = True
                                        break
                            
                            # If no segment information found, use the provided segment or mark as unknown
                            if not segment_found:
                                if segment and segment != "all":
                                    metadata["segment"] = segment
                                    segment_counts[segment] += 1
                                else:
                                    metadata["segment"] = "unknown"
                                    segment_counts["unknown"] += 1
                    
                    # Only include records matching the requested segment if specified
                    if segment and segment != "all" and metadata["segment"] != segment:
                        continue
                    
                    # Add to sequences list
                    sequences.append(record)
                    
                    # Add to metadata lines
                    metadata_line = "\t".join([
                        metadata["strain"],
                        metadata["virus"],
                        metadata["accession"],
                        metadata["date"],
                        metadata["country"],
                        metadata["division"],
                        metadata["location"],
                        metadata["host"],
                        metadata["segment"],
                        str(metadata["length"])
                    ])
                    metadata_lines.append(metadata_line)
                
                # Sleep between batches to avoid overwhelming NCBI
                if end < len(id_list):
                    time.sleep(1)
            
            # Log results
            log_with_context(logger, "INFO", f"Successfully fetched {len(sequences)} sequences", 
                           {"segment_counts": segment_counts})
            
            return {
                "sequences": sequences,
                "metadata_lines": metadata_lines,
                "stats": {
                    "total": len(sequences),
                    "segments": segment_counts,
                    "last_id": last_id
                }
            }
            
        except Exception as e:
            log_with_context(logger, "ERROR", f"Error fetching sequences: {e}")
            return {"sequences": [], "metadata_lines": ["strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength"], "stats": {"total": 0, "segments": {}}}
    
    def fetch_sequences(self, query: str, max_sequences: int = 500, from_date: Optional[str] = None, segment: Optional[str] = None) -> Dict:
        """
        Search and fetch sequences matching the query
        
        Args:
            query: Search term for NCBI
            max_sequences: Maximum number of sequences to return
            from_date: Only return results from this date forward (YYYY-MM-DD)
            segment: Filter by specific segment (L, M, S, or all)
            
        Returns:
            Dictionary with sequences, metadata lines, and stats
        """
        # Search for sequence IDs
        id_list, web_env, query_key = self._search_sequences(query, max_sequences, from_date)
        
        # If no results, return empty dict
        if not id_list:
            log_with_context(logger, "INFO", "No sequences found matching search criteria")
            return {"sequences": [], "metadata_lines": ["strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength"], "stats": {"total": 0, "segments": {}}}
        
        # Fetch sequences
        return self._fetch_sequences(id_list, web_env, query_key, segment)

class DataArchiver:
    """Manages versioned archiving of sequence data"""
    
    def __init__(self, archive_dir: str):
        self.archive_dir = archive_dir
        log_with_context(logger, "INFO", f"Initialized data archiver with directory: {archive_dir}")
    
    def archive_data(self, sequences_file: str, metadata_file: str, pathogen: str) -> Dict[str, str]:
        """
        Archive sequence and metadata files with versioning
        
        Args:
            sequences_file: Path to sequences FASTA file
            metadata_file: Path to metadata TSV file
            pathogen: Pathogen name (used for directory naming)
            
        Returns:
            Dictionary with paths to archived files
        """
        try:
            # Create timestamp for versioning
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            
            # Create archive paths
            archive_base = os.path.join(self.archive_dir, "archive", pathogen, timestamp)
            os.makedirs(archive_base, exist_ok=True)
            
            sequences_archive = os.path.join(archive_base, f"{pathogen}_sequences.fasta")
            metadata_archive = os.path.join(archive_base, f"{pathogen}_metadata.tsv")
            
            # Copy files to archive
            shutil.copy2(sequences_file, sequences_archive)
            shutil.copy2(metadata_file, metadata_archive)
            
            log_with_context(logger, "INFO", f"Archived data to {archive_base}", 
                           {"sequences": sequences_archive, "metadata": metadata_archive})
            
            return {"sequences": sequences_archive, "metadata": metadata_archive}
            
        except Exception as e:
            log_with_context(logger, "ERROR", f"Error archiving data: {e}", 
                           {"sequences_file": sequences_file, "metadata_file": metadata_file})
            return {}

def parse_args():
    parser = argparse.ArgumentParser(description="Download sequences and metadata for pathogen genomic analysis")
    
    parser.add_argument("--search-term", required=True, help="Search term for NCBI Entrez")
    parser.add_argument("--max-sequences", type=int, default=500, help="Maximum number of sequences to download")
    parser.add_argument("--output-sequences", required=True, help="Output path for sequences FASTA file")
    parser.add_argument("--output-metadata", required=True, help="Output path for metadata TSV file")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--segment", default="", help="Filter by specific segment (L, M, S, or all)")
    parser.add_argument("--incremental", action="store_true", help="Only download sequences added since the last run")
    parser.add_argument("--archive", action="store_true", help="Archive downloaded data with timestamps")
    parser.add_argument("--archive-dir", default="data", help="Directory for archived data")
    parser.add_argument("--tracking-file", default="config/data_tracking.json", help="Path to data tracking file")
    parser.add_argument("--pathogen", default="pathogen", help="Pathogen name for archiving")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
                        help="Logging level")
    
    return parser.parse_args()

def main():
    start_time = time.time()
    args = parse_args()
    
    # Update logger level based on arguments
    logger.setLevel(logging.getLevelName(args.log_level))
    
    try:
        # Initialize data tracker
        tracker = DataTracker(args.tracking_file)
        
        # Set up data source
        ncbi_source = NCBIDataSource(args.email)
        
        # Determine from_date for incremental fetching
        from_date = None
        if args.incremental:
            from_date = tracker.get_last_update_date()
            if from_date:
                log_with_context(logger, "INFO", f"Performing incremental update from {from_date}")
            else:
                log_with_context(logger, "WARNING", "No previous update date found, performing full download")
        
        # Fetch sequences
        fetch_results = ncbi_source.fetch_sequences(
            args.search_term,
            args.max_sequences,
            from_date=from_date,
            segment=args.segment
        )
        
        if fetch_results["stats"]["total"] == 0:
            log_with_context(logger, "INFO", "No sequences found or no new sequences since last update")
            log_execution_stats(logger, start_time, {"operation": "download", "sequences": 0, "status": "no_data"})
            return
        
        # Create output directories
        os.makedirs(os.path.dirname(args.output_sequences), exist_ok=True)
        os.makedirs(os.path.dirname(args.output_metadata), exist_ok=True)
        
        # Write sequences to FASTA
        with open(args.output_sequences, "w") as f:
            SeqIO.write(fetch_results["sequences"], f, "fasta")
        
        # Write metadata to TSV
        with open(args.output_metadata, "w") as f:
            f.write("\n".join(fetch_results["metadata_lines"]))
        
        # Archive data if enabled
        archive_paths = {}
        if args.archive:
            archiver = DataArchiver(args.archive_dir)
            archive_paths = archiver.archive_data(
                args.output_sequences,
                args.output_metadata,
                args.pathogen
            )
        
        # Update tracking data
        tracker.update({
            "total": fetch_results["stats"]["total"],
            "segments": fetch_results["stats"]["segments"],
            "last_id": fetch_results["stats"].get("last_id")
        })
        
        log_with_context(logger, "INFO", "Data download completed successfully", {
            "sequences": fetch_results["stats"]["total"],
            "segments": fetch_results["stats"]["segments"],
            "files": {
                "sequences": args.output_sequences,
                "metadata": args.output_metadata
            }
        })
        
        # Log execution statistics
        log_execution_stats(logger, start_time, {
            "operation": "download",
            "sequences": fetch_results["stats"]["total"],
            "segments": fetch_results["stats"]["segments"]
        })
        
    except Exception as e:
        log_with_context(logger, "CRITICAL", f"Unhandled exception: {e}")
        log_execution_stats(logger, start_time, {"operation": "download", "error": str(e)}, status="failed")
        sys.exit(1)

if __name__ == "__main__":
    main()