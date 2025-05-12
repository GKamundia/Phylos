#!/usr/bin/env python3
"""
Enhanced sequence download utility for pathogen genomic surveillance
Supports configurable data sources, versioned archiving, and incremental updates
"""

import os
import sys
import json
import snakemake
import yaml
import argparse
import datetime
import time
import logging
from pathlib import Path
import shutil
from typing import Dict, List, Optional, Union
from Bio import Entrez
from Bio import SeqIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

class DataTracker:
    """Tracks downloaded data for incremental updates"""
    
    def __init__(self, tracking_file_path):
        self.tracking_file_path = tracking_file_path
        self._load_tracking_data()
    
    def _load_tracking_data(self):
        """Load tracking data from JSON file"""
        if os.path.exists(self.tracking_file_path):
            with open(self.tracking_file_path, 'r') as f:
                try:
                    self.tracking_data = json.load(f)
                    logger.info(f"Loaded tracking data from {self.tracking_file_path}")
                except json.JSONDecodeError:
                    logger.error(f"Error parsing tracking file {self.tracking_file_path}")
                    self._initialize_tracking_data()
        else:
            logger.warning(f"Tracking file {self.tracking_file_path} not found. Creating new tracking data.")
            self._initialize_tracking_data()
    
    def _initialize_tracking_data(self):
        """Initialize tracking data structure"""
        today = datetime.datetime.now().strftime("%Y-%m-%d")
        self.tracking_data = {
            "last_update": today,
            "sequences": {
                "total_count": 0,
                "last_sequence_id": None
            },
            "segments": {
                "L": {"count": 0, "last_update": None},
                "M": {"count": 0, "last_update": None},
                "S": {"count": 0, "last_update": None}
            }
        }
    
    def save(self):
        """Save tracking data to JSON file"""
        with open(self.tracking_file_path, 'w') as f:
            json.dump(self.tracking_data, f, indent=2)
            logger.info(f"Updated tracking data saved to {self.tracking_file_path}")
    
    def get_last_update_date(self) -> str:
        """Get the last update date"""
        return self.tracking_data.get("last_update", "")
    
    def get_segment_count(self, segment: str) -> int:
        """Get sequence count for a specific segment"""
        if segment in self.tracking_data["segments"]:
            return self.tracking_data["segments"][segment].get("count", 0)
        return 0
    
    def get_total_count(self) -> int:
        """Get total sequence count"""
        return self.tracking_data["sequences"].get("total_count", 0)
    
    def update_after_download(self, download_stats: Dict):
        """Update tracking data after a download session"""
        today = datetime.datetime.now().strftime("%Y-%m-%d")
        self.tracking_data["last_update"] = today
        self.tracking_data["sequences"]["total_count"] += download_stats["total_added"]
        
        # Update segment-specific stats
        for segment, count in download_stats["segments"].items():
            if segment in self.tracking_data["segments"]:
                self.tracking_data["segments"][segment]["count"] += count
                if count > 0:
                    self.tracking_data["segments"][segment]["last_update"] = today
            else:
                # Handle new segment (flexibility for other pathogens)
                self.tracking_data["segments"][segment] = {
                    "count": count,
                    "last_update": today if count > 0 else None
                }
        
        # Update last sequence ID if provided
        if download_stats.get("last_sequence_id"):
            self.tracking_data["sequences"]["last_sequence_id"] = download_stats["last_sequence_id"]
        
        self.save()

class NCBIDataSource:
    """Handles sequence fetching from NCBI GenBank"""
    
    def __init__(self, email: str, retry_attempts: int = 3, retry_delay: int = 5):
        self.email = email
        self.retry_attempts = retry_attempts
        self.retry_delay = retry_delay
        Entrez.email = email
        logger.info(f"Initialized NCBI data source with email: {email}")
    
    def fetch_sequences(self, search_term: str, max_sequences: int, 
                        from_date: Optional[str] = None, segment: Optional[str] = None) -> Dict:
        """
        Fetch sequences from NCBI based on search criteria
        
        Args:
            search_term: NCBI search query
            max_sequences: Maximum number of sequences to retrieve
            from_date: Optional date to fetch sequences from (YYYY/MM/DD)
            segment: Optional segment filter
            
        Returns:
            Dictionary with sequences and metadata
        """
        # Add date range to search term if specified
        if from_date:
            # Convert to NCBI date format if needed
            try:
                date_obj = datetime.datetime.strptime(from_date, "%Y-%m-%d")
                ncbi_date = date_obj.strftime("%Y/%m/%d")
                search_term = f"{search_term} AND {ncbi_date}[PDAT] : 3000[PDAT]"
                logger.info(f"Searching for sequences from {from_date} onwards")
            except ValueError:
                logger.warning(f"Invalid date format: {from_date}. Expected YYYY-MM-DD.")
        
        # Search for sequences
        logger.info(f"Searching NCBI for: {search_term}")
        for attempt in range(1, self.retry_attempts + 1):
            try:
                search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=max_sequences)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                break
            except Exception as e:
                if attempt < self.retry_attempts:
                    logger.warning(f"Search attempt {attempt} failed: {e}. Retrying in {self.retry_delay} seconds...")
                    time.sleep(self.retry_delay)
                else:
                    logger.error(f"Failed to search NCBI after {self.retry_attempts} attempts: {e}")
                    raise
        
        id_list = search_results["IdList"]
        if not id_list:
            logger.warning("No sequences found matching search criteria")
            return {"sequences": [], "metadata_lines": [], "stats": {"total": 0}}
        
        logger.info(f"Found {len(id_list)} sequences in NCBI")
        
        # Download sequences in batches
        batch_size = 100
        sequences = []
        metadata_lines = ["strain\tvirus\taccession\tdate\tcountry\tdivision\tlocation\thost\tsegment\tlength"]
        segment_counts = {"L": 0, "M": 0, "S": 0, "unknown": 0}
        last_id = None
        
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i+batch_size]
            logger.info(f"Downloading batch {i//batch_size + 1} of {(len(id_list)-1)//batch_size + 1}")
            
            for attempt in range(1, self.retry_attempts + 1):
                try:
                    fetch_handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="text")
                    break
                except Exception as e:
                    if attempt < self.retry_attempts:
                        logger.warning(f"Fetch attempt {attempt} failed: {e}. Retrying in {self.retry_delay} seconds...")
                        time.sleep(self.retry_delay)
                    else:
                        logger.error(f"Failed to fetch records after {self.retry_attempts} attempts: {e}")
                        raise
            
            for record in SeqIO.parse(fetch_handle, "genbank"):
                # Extract metadata
                accession = record.id
                last_id = accession  # Keep track of last ID for resuming
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
                if segment and segment_info and segment.lower() != segment_info.lower():
                    continue
                
                # Track segment counts
                if segment_info:
                    if segment_info.upper() in segment_counts:
                        segment_counts[segment_info.upper()] += 1
                    else:
                        segment_counts["unknown"] += 1
                else:
                    segment_counts["unknown"] += 1
                
                # Add to collections
                sequences.append(record)
                metadata_lines.append(f"{strain}\tPathogen\t{accession}\t{date}\t{country}\t{division}\t\t{host}\t{segment_info}\t{length}")
            
            fetch_handle.close()
        
        return {
            "sequences": sequences,
            "metadata_lines": metadata_lines,
            "stats": {
                "total": len(sequences),
                "segments": segment_counts,
                "last_id": last_id
            }
        }

class DataArchiver:
    """Manages versioned archiving of sequence data"""
    
    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self.today = datetime.datetime.now().strftime("%Y-%m-%d")
        
    def get_archive_paths(self, pathogen: str) -> Dict[str, Path]:
        """Generate archive paths for the current date"""
        archive_date_dir = self.base_dir / "archive" / self.today
        
        paths = {
            "sequences_dir": archive_date_dir / "sequences",
            "metadata_dir": archive_date_dir / "metadata",
            "sequences_file": archive_date_dir / "sequences" / f"{pathogen}_sequences_{self.today}.fasta",
            "metadata_file": archive_date_dir / "metadata" / f"{pathogen}_metadata_{self.today}.tsv"
        }
        
        # Ensure directories exist
        paths["sequences_dir"].mkdir(parents=True, exist_ok=True)
        paths["metadata_dir"].mkdir(parents=True, exist_ok=True)
        
        return paths
        
    def archive_data(self, sequences_file: str, metadata_file: str, pathogen: str) -> Dict[str, str]:
        """
        Archive sequence and metadata files
        
        Args:
            sequences_file: Path to original sequences FASTA
            metadata_file: Path to original metadata TSV
            pathogen: Pathogen name for file naming
            
        Returns:
            Dictionary with archive file paths
        """
        paths = self.get_archive_paths(pathogen)
        
        # Copy files to archive
        logger.info(f"Archiving sequences to {paths['sequences_file']}")
        shutil.copy2(sequences_file, paths['sequences_file'])
        
        logger.info(f"Archiving metadata to {paths['metadata_file']}")
        shutil.copy2(metadata_file, paths['metadata_file'])
        
        return {
            "sequences": str(paths['sequences_file']),
            "metadata": str(paths['metadata_file'])
        }

def parse_args():
    parser = argparse.ArgumentParser(
        description="Enhanced pathogen sequence downloader with archiving and incremental updates")
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
    parser.add_argument('--tracking-file', default="config/data_tracking.json",
                        help="Path to data tracking JSON file")
    parser.add_argument('--archive', action='store_true',
                        help="Enable versioned archiving of downloaded data")
    parser.add_argument('--archive-dir', default="data",
                        help="Base directory for archived files")
    parser.add_argument('--pathogen', default="pathogen",
                        help="Pathogen name for archive file naming")
    parser.add_argument('--incremental', action='store_true',
                        help="Enable incremental updates (fetch only new sequences)")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Initialize data tracker
    tracker = DataTracker(args.tracking_file)
    
    # Set up data source
    ncbi_source = NCBIDataSource(args.email)
    
    # Determine from_date for incremental fetching
    from_date = None
    if args.incremental:
        from_date = tracker.get_last_update_date()
        if from_date:
            logger.info(f"Performing incremental update from {from_date}")
        else:
            logger.warning("No previous update date found, performing full download")
    
    # Fetch sequences
    fetch_results = ncbi_source.fetch_sequences(
        args.search_term,
        args.max_sequences,
        from_date=from_date,
        segment=args.segment
    )
    
    if fetch_results["stats"]["total"] == 0:
        logger.info("No sequences found or no new sequences since last update")
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
    if args.archive:
        archiver = DataArchiver(args.archive_dir)
        archive_paths = archiver.archive_data(
            args.output_sequences,
            args.output_metadata,
            args.pathogen
        )
        logger.info(f"Data archived to {os.path.dirname(archive_paths['sequences'])}")
    
    # Update tracking data
    download_stats = {
        "total_added": fetch_results["stats"]["total"],
        "segments": fetch_results["stats"]["segments"],
        "last_sequence_id": fetch_results["stats"].get("last_id")
    }
    tracker.update_after_download(download_stats)
    
    logger.info(f"Downloaded {fetch_results['stats']['total']} sequences")
    for segment, count in fetch_results["stats"]["segments"].items():
        if count > 0:
            logger.info(f"  - Segment {segment}: {count} sequences")
    
    logger.info(f"Sequences saved to: {args.output_sequences}")
    logger.info(f"Metadata saved to: {args.output_metadata}")

if __name__ == "__main__":
    main()

# scripts/run_download.py
import os
import subprocess
import sys

# Get parameters from Snakemake
output_sequences = snakemake.output.sequences
output_metadata = snakemake.output.metadata
email = snakemake.params.email
search_term = snakemake.params.search_term
max_sequences = snakemake.params.max_sequences
segment = snakemake.params.segment
tracking_file = snakemake.params.tracking_file
archive = snakemake.params.archive
pathogen = snakemake.params.pathogen
incremental = snakemake.params.incremental
log_file = snakemake.log[0]

# Build command
cmd = [
    "python", "scripts/download_sequences.py",
    "--search-term", search_term,
    "--max-sequences", str(max_sequences),
    "--output-sequences", output_sequences,
    "--output-metadata", output_metadata,
    "--email", email,
    "--tracking-file", tracking_file,
    "--archive-dir", "data",
    "--pathogen", pathogen
]

# Add conditional arguments
if segment:
    cmd.extend(["--segment", segment])
if archive:
    cmd.append("--archive")
if incremental:
    cmd.append("--incremental")

# Execute command
with open(log_file, "w") as log:
    result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    
sys.exit(result.returncode)