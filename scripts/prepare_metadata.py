#!/usr/bin/env python3
"""
Enhanced metadata preparation script for Rift Valley Fever Nextstrain analysis.
- Validates metadata against formal JSON schema
- Standardizes dates, locations, host names
- Infers missing metadata where possible
- Reports quality metrics and flags problematic records
"""

import os
import sys
import json
import pandas as pd
import argparse
import re
import logging
from datetime import datetime
from pathlib import Path
import jsonschema
from jsonschema import validate, ValidationError

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Host name standardization mapping
HOST_MAPPING = {
    "human": "Homo sapiens",
    "humans": "Homo sapiens",
    "homo sapiens": "Homo sapiens",
    "h. sapiens": "Homo sapiens",
    "sheep": "Ovis aries",
    "domestic sheep": "Ovis aries",
    "ovis aries": "Ovis aries",
    "cow": "Bos taurus",
    "cows": "Bos taurus",
    "cattle": "Bos taurus",
    "bovine": "Bos taurus",
    "bos taurus": "Bos taurus",
    "goat": "Capra hircus",
    "goats": "Capra hircus",
    "capra": "Capra hircus",
    "capra hircus": "Capra hircus",
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Enhanced metadata preparation for RVF Nextstrain analysis")
    parser.add_argument('input', help="Path to input metadata TSV file")
    parser.add_argument('output', help="Path to output metadata TSV file")
    parser.add_argument('--schema', default="config/metadata_schema.json", 
                        help="Path to JSON schema file")
    parser.add_argument('--lat-longs', default="config/lat_longs.tsv", 
                        help="Path to lat_longs.tsv file")
    parser.add_argument('--report', default=None,
                        help="Path to write validation report (optional)")
    parser.add_argument('--strict', action='store_true',
                        help="Fail if any records don't pass validation")
    return parser.parse_args()

def load_schema(schema_file):
    """Load and parse the JSON schema file"""
    try:
        with open(schema_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load schema file: {e}")
        sys.exit(1)

def load_lat_longs(lat_longs_file):
    """Load latitude/longitude reference data"""
    try:
        df = pd.read_csv(lat_longs_file, sep='\t', comment='#')
        # Convert to a nested dictionary for easier lookups
        result = {}
        for _, row in df.iterrows():
            location_type = row[0]  # country, division, etc.
            location_name = row[1]
            lat = row[2]
            long = row[3]
            
            if location_type not in result:
                result[location_type] = {}
            
            result[location_type][location_name] = (lat, long)
        
        return result
    except Exception as e:
        logger.error(f"Failed to load lat_longs file: {e}")
        return {}

def standardize_dates(date_str):
    """Convert various date formats to YYYY-MM-DD, YYYY-MM, or YYYY"""
    if not date_str or pd.isna(date_str):
        return ""
    
    # Try common date formats
    date_formats = [
        "%Y-%m-%d", "%Y/%m/%d",  # YYYY-MM-DD, YYYY/MM/DD
        "%d-%m-%Y", "%d/%m/%Y",  # DD-MM-YYYY, DD/MM/YYYY
        "%m-%d-%Y", "%m/%d/%Y",  # MM-DD-YYYY, MM/DD/YYYY
        "%Y-%m", "%Y/%m",        # YYYY-MM, YYYY/MM
        "%b %Y", "%B %Y",        # MMM YYYY (e.g. Jan 2020)
        "%d %b %Y", "%d %B %Y",  # DD MMM YYYY (e.g. 15 Jan 2020)
        "%b %d %Y", "%B %d %Y",  # MMM DD YYYY (e.g. Jan 15 2020)
        "%Y"                     # YYYY
    ]
    
    # Try to parse the date
    for fmt in date_formats:
        try:
            date_obj = datetime.strptime(str(date_str).strip(), fmt)
            # Format based on precision of the input
            if fmt == "%Y":
                return f"{date_obj.year}"
            elif fmt in ["%Y-%m", "%Y/%m", "%b %Y", "%B %Y"]:
                return f"{date_obj.year}-{date_obj.month:02d}"
            else:
                return date_obj.strftime("%Y-%m-%d")
        except ValueError:
            continue
    
    # Extract year from string (last resort)
    year_match = re.search(r'(\d{4})', str(date_str))
    if year_match:
        return year_match.group(1)
    
    return ""

def standardize_host(host_str):
    """Normalize host names to standard scientific names"""
    if not host_str or pd.isna(host_str):
        return ""
    
    host_lower = str(host_str).lower().strip()
    
    # Check direct mapping
    if host_lower in HOST_MAPPING:
        return HOST_MAPPING[host_lower]
    
    # Keep original if no mapping found
    return host_str

def infer_division(row, lat_longs):
    """Attempt to infer division from location or other fields"""
    # This is a placeholder for more complex inference logic
    # For now, we just check if there's a division in the lat_longs file
    
    if pd.isna(row.division) or not row.division:
        if not pd.isna(row.country) and row.country:
            # Check if there are divisions for this country in lat_longs
            if 'division' in lat_longs:
                divisions_for_country = [div for div, (lat, long) in lat_longs['division'].items() 
                                        if div.startswith(f"{row.country}:")]
                if len(divisions_for_country) == 1:
                    # If there's only one division for this country, it's likely the capital or main region
                    return divisions_for_country[0].split(':')[1]
    
    return row.division

def get_coordinates(location_type, location_name, lat_longs):
    """Get coordinates for a location from lat_longs dictionary"""
    if location_type in lat_longs and location_name in lat_longs[location_type]:
        return lat_longs[location_type][location_name]
    return None

def validate_record(record, schema):
    """Validate a single metadata record against the schema"""
    try:
        validate(instance=record, schema=schema)
        return True, None
    except ValidationError as e:
        return False, str(e)

def validate_metadata(metadata, schema):
    """Validate all metadata records against the schema"""
    validation_results = []
    
    for idx, row in metadata.iterrows():
        record = row.to_dict()
        valid, error = validate_record(record, schema)
        
        if not valid:
            validation_results.append({
                'row': idx,
                'strain': row.get('strain', 'UNKNOWN'),
                'valid': False,
                'error': error
            })
    
    return validation_results

def generate_report(metadata, validation_results, report_file):
    """Generate a validation and quality report"""
    total_records = len(metadata)
    invalid_records = len(validation_results)
    
    # Calculate completion percentages for key fields
    completion = {}
    for field in ['strain', 'date', 'country', 'division', 'host', 'segment']:
        if field in metadata.columns:
            non_empty = metadata[field].notna() & (metadata[field] != '')
            completion[field] = f"{non_empty.sum() / total_records * 100:.1f}%"
    
    # Create report content
    report = {
        'summary': {
            'total_records': total_records,
            'valid_records': total_records - invalid_records,
            'invalid_records': invalid_records,
            'completion_rate': completion
        },
        'invalid_records': validation_results
    }
    
    # Write report to file
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    return report

def main():
    args = parse_args()
    
    # Load schema
    schema = load_schema(args.schema)
    logger.info(f"Loaded schema from {args.schema}")
    
    # Load lat_longs reference data
    lat_longs = load_lat_longs(args.lat_longs)
    logger.info(f"Loaded geographic coordinates from {args.lat_longs}")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Load metadata
    logger.info(f"Loading metadata from {args.input}")
    try:
        metadata = pd.read_csv(args.input, sep='\t', dtype=str)
    except Exception as e:
        logger.error(f"Error loading metadata: {e}")
        return 1
    
    logger.info(f"Loaded {len(metadata)} records")
    
    # Clean and standardize metadata
    logger.info("Standardizing metadata...")
    
    # Standardize dates
    if 'date' in metadata.columns:
        metadata['date'] = metadata['date'].apply(standardize_dates)
        logger.info(f"Standardized dates: {metadata['date'].nunique()} unique values")
    
    # Standardize hosts
    if 'host' in metadata.columns:
        metadata['host'] = metadata['host'].apply(standardize_host)
        logger.info(f"Standardized hosts: {metadata['host'].nunique()} unique hosts")
    
    # Ensure strain names are unique
    if 'strain' in metadata.columns:
        # If duplicates exist, append a suffix
        duplicated = metadata['strain'].duplicated()
        if any(duplicated):
            logger.warning(f"Found {sum(duplicated)} duplicate strain names")
            dup_strains = metadata.loc[duplicated, 'strain'].unique()
            for strain in dup_strains:
                mask = metadata['strain'] == strain
                indices = metadata.index[mask]
                for i, idx in enumerate(indices[1:], start=1):
                    metadata.loc[idx, 'strain'] = f"{strain}_{i}"
            logger.info("Added suffix to duplicate strain names")
    
    # Infer missing metadata where possible
    if 'division' in metadata.columns and 'country' in metadata.columns:
        # Apply inference function to each row
        mask = metadata['division'].isna() | (metadata['division'] == '')
        if any(mask):
            logger.info(f"Attempting to infer division for {sum(mask)} records")
            metadata.loc[mask, 'division'] = metadata[mask].apply(
                lambda row: infer_division(row, lat_longs), axis=1
            )
    
    # Check for required fields
    required_fields = schema.get('required', [])
    for field in required_fields:
        if field not in metadata.columns:
            logger.error(f"Required field '{field}' is missing from metadata")
            if args.strict:
                return 1
        else:
            missing = metadata[field].isna() | (metadata[field] == '')
            if any(missing):
                logger.warning(f"{sum(missing)} records missing required field '{field}'")
    
    # Validate against schema
    logger.info("Validating metadata against schema...")
    validation_results = validate_metadata(metadata, schema)
    
    if validation_results:
        logger.warning(f"Found {len(validation_results)} records with validation errors")
        if args.strict:
            logger.error("Validation failed in strict mode")
            return 1
    else:
        logger.info("All records passed schema validation")
    
    # Generate report if requested
    if args.report:
        report_file = args.report
        report = generate_report(metadata, validation_results, report_file)
        logger.info(f"Validation report written to {report_file}")
    
    # Save cleaned metadata
    logger.info(f"Saving cleaned metadata to {args.output}")
    metadata.to_csv(args.output, sep='\t', index=False)
    logger.info("Metadata preparation completed successfully")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())