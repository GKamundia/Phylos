#!/usr/bin/env python3
"""
Fix segment identification for RVF virus sequences
Specifically handles glycoprotein and other M segment indicators
"""

import pandas as pd
import re
from pathlib import Path
import logging
import os

# Check if running under Snakemake
try:
    import snakemake
    SNAKEMAKE_MODE = True
    
    # Get input parameters from Snakemake
    sequences_file = snakemake.input.sequences
    metadata_file = snakemake.input.metadata
    output_metadata_file = snakemake.output.metadata
    log_file = snakemake.log[0]
    
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)
    
except (ImportError, AttributeError):
    SNAKEMAKE_MODE = False
    metadata_file = "results/filtered/rvf_metadata.tsv"
    output_metadata_file = metadata_file
    
    # Set up basic logging for standalone mode
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)

def identify_rvf_segment(title, length=None):
    """
    Enhanced segment identification for RVF virus
    
    Args:
        title (str): GenBank title/description
        length (int): Sequence length (optional, for additional validation)
    
    Returns:
        str: Segment (L, M, S) or empty string if unclear
    """
    if not title or pd.isna(title):
        return ""
    
    title_lower = str(title).lower()
    
    # Direct segment mentions (highest priority)
    if re.search(r'\b(segment\s*s|s\s*segment)\b', title_lower):
        return "S"
    elif re.search(r'\b(segment\s*m|m\s*segment)\b', title_lower):
        return "M"  
    elif re.search(r'\b(segment\s*l|l\s*segment)\b', title_lower):
        return "L"
    
    # M segment indicators (glycoproteins)
    m_indicators = [
        r'\bglycoprote?in\b',  # glycoprotein or glycoprotein
        r'\bg[nc12]\b',        # Gn, Gc, G1, G2
        r'\benvelope\s+protein\b',
        r'\bsurface\s+protein\b',
        r'\bmembrane\s+protein\b',
        r'\bglycosylated\s+protein\b'
    ]
    
    for pattern in m_indicators:
        if re.search(pattern, title_lower):
            return "M"
    
    # L segment indicators (polymerase)
    l_indicators = [
        r'\bpolymerase\b',
        r'\brna\s+polymerase\b',
        r'\bl\s+protein\b',
        r'\blarge\s+protein\b',
        r'\breplicase\b'
    ]
    
    for pattern in l_indicators:
        if re.search(pattern, title_lower):
            return "L"
    
    # S segment indicators (nucleocapsid, NSs)
    s_indicators = [
        r'\bnucleocapsid\b',
        r'\bn\s+protein\b',
        r'\bnss?\s+protein\b',
        r'\bsmall\s+protein\b',
        r'\bstructural\s+protein\b'
    ]
    
    for pattern in s_indicators:
        if re.search(pattern, title_lower):
            return "S"
    


def fix_segment_identification():
    """Fix segment identification in the filtered metadata"""
    
    # Load the filtered metadata
    logger.info(f"Loading metadata from {metadata_file}")
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    logger.info(f"Loaded {len(metadata)} records from {metadata_file}")
      # Count current segment distribution
    current_segments = metadata['Segment'].value_counts(dropna=False)
    logger.info(f"Current segment distribution:")
    for segment, count in current_segments.items():
        logger.info(f"  {segment}: {count}")
    
    # Create a standardized lowercase 'segment' column from the uppercase 'Segment' column
    if 'Segment' in metadata.columns:
        # Convert uppercase 'Segment' to lowercase 'segment'
        metadata['segment'] = metadata['Segment'].str.lower()
        logger.info("Created lowercase 'segment' column from 'Segment' column")
    elif 'segment' not in metadata.columns:
        # Create segment column if it doesn't exist
        metadata['segment'] = None
        logger.info("Created new 'segment' column")
    
    # Find records with missing segments (check both uppercase and lowercase)
    missing_segments = metadata['segment'].isna()
    missing_count = missing_segments.sum()
    logger.info(f"Records with missing segments: {missing_count}")
    
    if missing_count > 0:
        logger.info("Attempting to fix missing segments...")
        
        # Fix missing segments
        fixed_count = 0
        for idx, row in metadata[missing_segments].iterrows():
            title = row.get('GenBank_Title', '')
            length = row.get('Length', None)
              # Try to identify segment
            segment = identify_rvf_segment(title, length)
            
            if segment:
                metadata.loc[idx, 'segment'] = segment.lower()  # Ensure lowercase
                fixed_count += 1
                logger.info(f"  Fixed {row['accession']}: {segment} (from: {title[:80]}...)")
        
        logger.info(f"Fixed {fixed_count} out of {missing_count} missing segments")
          # Show new distribution (using lowercase segment column)
        new_segments = metadata['segment'].value_counts(dropna=False)
        logger.info(f"New segment distribution:")
        for segment, count in new_segments.items():
            logger.info(f"  {segment}: {count}")
        
        # Save the corrected metadata
        if SNAKEMAKE_MODE:
            # In Snakemake mode, write to the designated output file
            os.makedirs(os.path.dirname(output_metadata_file), exist_ok=True)
            metadata.to_csv(output_metadata_file, sep='\t', index=False)
            logger.info(f"Updated metadata saved to: {output_metadata_file}")
        else:
            # In standalone mode, backup and overwrite original
            backup_file = metadata_file + ".backup"
            metadata.to_csv(backup_file, sep='\t', index=False)
            logger.info(f"Backup saved to: {backup_file}")
            
            metadata.to_csv(metadata_file, sep='\t', index=False)
            logger.info(f"Updated metadata saved to: {metadata_file}")
        
        return fixed_count
    else:
        logger.info("No missing segments to fix!")
        # Always create output file in Snakemake mode, even if no fixes needed
        if SNAKEMAKE_MODE:
            os.makedirs(os.path.dirname(output_metadata_file), exist_ok=True)
            metadata.to_csv(output_metadata_file, sep='\t', index=False)
            logger.info(f"Metadata copied to: {output_metadata_file}")
        return 0

if __name__ == "__main__":
    fixed_count = fix_segment_identification()
    if SNAKEMAKE_MODE:
        logger.info(f"✅ Segment identification fix complete: {fixed_count} records corrected")
    else:
        print(f"\n✅ Segment identification fix complete: {fixed_count} records corrected")
