#!/usr/bin/env python3
"""
Metadata preparation and standardization script for Nextstrain RVF analysis
"""

import logging
import os
import sys
import json
import time
import argparse
from datetime import datetime
from pathlib import Path
import pandas as pd
import re
from jsonschema import validate, ValidationError

# Configure basic logging to stderr for Snakemake to capture
logging.basicConfig(stream=sys.stderr, level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)

# Import our new logging utilities
sys.path.append(os.path.join(os.path.dirname(__file__), '.'))
from utils.log_utils import setup_logger, log_execution_stats, log_with_context

# Configure logger
logger = setup_logger(
    name="prepare_metadata", 
    log_file="prepare_metadata.log", 
    level="INFO",
    json_output=True
)

def parse_args():
    parser = argparse.ArgumentParser(description="Prepare and standardize metadata for Nextstrain analysis")
    
    parser.add_argument("input", help="Input metadata TSV file")
    parser.add_argument("output", help="Output metadata TSV file")
    parser.add_argument("--schema", help="JSON schema for metadata validation")
    parser.add_argument("--lat-longs", help="TSV file with latitude/longitude data")
    parser.add_argument("--strict", action="store_true", help="Enforce strict schema validation")
    parser.add_argument("--report", help="Output path for validation report JSON")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
                        help="Logging level")
    
    return parser.parse_args()

def load_schema(schema_file):
    """Load JSON schema file for metadata validation"""
    try:
        if not schema_file or not os.path.exists(schema_file):
            log_with_context(logger, "WARNING", f"Schema file not found: {schema_file}, skipping validation")
            return None
        
        with open(schema_file) as f:
            schema = json.load(f)
            
        log_with_context(logger, "INFO", f"Loaded schema from {schema_file}", 
                        {"required_fields": schema.get("required", [])})
        return schema
    
    except Exception as e:
        log_with_context(logger, "ERROR", f"Error loading schema: {e}", {"file": schema_file})
        return None

def load_lat_longs(lat_longs_file):
    """Load latitude/longitude reference data"""
    try:
        if not lat_longs_file or not os.path.exists(lat_longs_file):
            log_with_context(logger, "WARNING", f"Lat/long file not found: {lat_longs_file}, geographic coordinates will not be added")
            return {}
        
        # Read lat_longs file
        lat_longs_df = pd.read_csv(lat_longs_file, sep='\t')
        
        # Convert to dictionary for faster lookups
        lat_longs = {}
        for _, row in lat_longs_df.iterrows():
            location_type = row['type']
            location = row['location']
            lat = row['latitude']
            long = row['longitude']
            
            if location_type not in lat_longs:
                lat_longs[location_type] = {}
            
            lat_longs[location_type][location] = (lat, long)
        
        log_with_context(logger, "INFO", f"Loaded geographic coordinates from {lat_longs_file}", 
                        {"locations": sum(len(v) for v in lat_longs.values())})
        
        return lat_longs
    
    except Exception as e:
        log_with_context(logger, "ERROR", f"Error loading lat/longs: {e}", {"file": lat_longs_file})
        return {}

def standardize_date(date_str):
    """
    Convert various date formats to YYYY-MM-DD
    If only year or year-month is provided, keep as is
    """
    if not date_str or pd.isna(date_str):
        return ""
    
    date_str = str(date_str).strip()
    
    # Check if already in YYYY-MM-DD format
    if re.match(r'^\d{4}-\d{2}-\d{2}$', date_str):
        return date_str
    
    # Check if in YYYY-MM format
    if re.match(r'^\d{4}-\d{2}$', date_str):
        return date_str
    
    # Check if just year
    if re.match(r'^\d{4}$', date_str):
        return date_str
    
    # Try various formats
    formats = [
        '%Y-%m-%d', '%Y/%m/%d', '%d-%m-%Y', '%d/%m/%Y', 
        '%m-%d-%Y', '%m/%d/%Y', '%d-%b-%Y', '%d %b %Y', 
        '%b %d %Y', '%B %d %Y', '%d %B %Y'
    ]
    
    for fmt in formats:
        try:
            date_obj = datetime.strptime(date_str, fmt)
            return date_obj.strftime('%Y-%m-%d')
        except ValueError:
            continue
    
    # If all formats fail, log and return original
    log_with_context(logger, "WARNING", f"Could not standardize date: {date_str}")
    return date_str

def validate_metadata(metadata_df, schema):
    """
    Validate metadata against JSON schema
    Returns validation statistics and list of invalid records
    """
    if schema is None:
        return {"valid": True, "total_records": len(metadata_df), "invalid_records": 0, "errors": []}
    
    total_records = len(metadata_df)
    invalid_records = []
    errors = []
    
    # Check required fields
    required_fields = schema.get("required", [])
    for field in required_fields:
        if field not in metadata_df.columns:
            errors.append(f"Required field '{field}' missing from metadata")
    
    # If missing required fields, return early
    if errors:
        log_with_context(logger, "ERROR", f"Metadata missing required fields: {', '.join(errors)}")
        return {
            "valid": False, 
            "total_records": total_records, 
            "invalid_records": total_records, 
            "errors": errors
        }
    
    # Validate each record
    for idx, row in metadata_df.iterrows():
        record_errors = []
        
        # Check required fields have values
        for field in required_fields:
            if pd.isna(row[field]) or str(row[field]).strip() == "":
                record_errors.append(f"Missing required value for '{field}'")
        
        # Check field types and patterns if specified in schema
        properties = schema.get("properties", {})
        for field, config in properties.items():
            if field in row and not pd.isna(row[field]) and row[field] != "":
                # Check enum values
                if "enum" in config:
                    if str(row[field]) not in config["enum"] and not (None in config["enum"] and pd.isna(row[field])):
                        record_errors.append(f"Value '{row[field]}' for field '{field}' not in allowed values: {config['enum']}")
                
                # Check field pattern
                if "pattern" in config and isinstance(row[field], str):
                    pattern = config["pattern"]
                    if not re.match(pattern, row[field]):
                        record_errors.append(f"Value '{row[field]}' for field '{field}' does not match pattern: {pattern}")
        
        # Record invalid record
        if record_errors:
            invalid_records.append({
                "index": idx,
                "errors": record_errors,
                "record": row.to_dict()
            })
    
    # Log validation results
    validation_result = {
        "valid": len(invalid_records) == 0,
        "total_records": total_records,
        "invalid_records": len(invalid_records),
        "errors": errors,
        "invalid_record_details": invalid_records[:10]  # Limit to first 10 for brevity
    }
    
    if len(invalid_records) > 0:
        log_with_context(logger, "WARNING", f"Found {len(invalid_records)} invalid records out of {total_records}", 
                       {"invalid_percentage": f"{(len(invalid_records)/total_records)*100:.1f}%"})
    else:
        log_with_context(logger, "INFO", f"All {total_records} records passed validation")
    
    return validation_result

def add_geographic_coordinates(metadata_df, lat_longs):
    """
    Add latitude and longitude to metadata based on location information
    """
    if not lat_longs:
        return metadata_df
    
    # Add lat/long columns if not present
    for col in ['latitude', 'longitude']:
        if col not in metadata_df.columns:
            metadata_df[col] = ""
    
    # Fill in coordinates
    records_with_coords = 0
    records_missing_coords = 0
    
    for idx, row in metadata_df.iterrows():
        coords_found = False
        
        # Try country + division + location
        if all(x in row and not pd.isna(row[x]) and row[x] != "" for x in ['country', 'division', 'location']):
            location_key = f"{row['country']}:{row['division']}:{row['location']}"
            if 'location' in lat_longs and location_key in lat_longs['location']:
                lat, long = lat_longs['location'][location_key]
                metadata_df.at[idx, 'latitude'] = lat
                metadata_df.at[idx, 'longitude'] = long
                coords_found = True
        
        # Try country + division
        if not coords_found and all(x in row and not pd.isna(row[x]) and row[x] != "" for x in ['country', 'division']):
            division_key = f"{row['country']}:{row['division']}"
            if 'division' in lat_longs and division_key in lat_longs['division']:
                lat, long = lat_longs['division'][division_key]
                metadata_df.at[idx, 'latitude'] = lat
                metadata_df.at[idx, 'longitude'] = long
                coords_found = True
        
        # Try just country
        if not coords_found and 'country' in row and not pd.isna(row['country']) and row['country'] != "":
            if 'country' in lat_longs and row['country'] in lat_longs['country']:
                lat, long = lat_longs['country'][row['country']]
                metadata_df.at[idx, 'latitude'] = lat
                metadata_df.at[idx, 'longitude'] = long
                coords_found = True
        
        if coords_found:
            records_with_coords += 1
        else:
            records_missing_coords += 1
    
    log_with_context(logger, "INFO", f"Added coordinates to {records_with_coords} records", 
                   {"missing_coords": records_missing_coords})
    
    return metadata_df

def main():
    start_time = time.time()
    args = parse_args()
    
    # Ensure the logger level from args is applied if basicConfig set a default
    if hasattr(args, 'log_level'):
         logger.setLevel(logging.getLevelName(args.log_level))
    
    try:
        # Load schema
        schema = load_schema(args.schema)
        
        # Load lat_longs reference data
        lat_longs = load_lat_longs(args.lat_longs)
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        
        # Load metadata
        log_with_context(logger, "INFO", f"Loading metadata from {args.input}")
        try:
            metadata = pd.read_csv(args.input, sep='\t', dtype=str)
        except Exception as e:
            log_with_context(logger, "ERROR", f"Error loading metadata: {e}")
            return 1
        
        log_with_context(logger, "INFO", f"Loaded {len(metadata)} records")
        
        # Clean and standardize metadata
        log_with_context(logger, "INFO", "Standardizing metadata...")
        
        # Standardize dates
        if 'date' in metadata.columns:
            metadata['date'] = metadata['date'].apply(standardize_date)
            date_stats = metadata['date'].str.len().value_counts()
            log_with_context(logger, "INFO", "Date standardization complete", {
                "complete_dates": int(date_stats.get(10, 0)),  # YYYY-MM-DD (length 10)
                "partial_dates": {
                    "year_month": int(date_stats.get(7, 0)),  # YYYY-MM (length 7)
                    "year": int(date_stats.get(4, 0))  # YYYY (length 4)
                },
                "empty_dates": int(metadata['date'].isna().sum())
            })
        
        # Add geographic coordinates
        metadata = add_geographic_coordinates(metadata, lat_longs)
        
        # Convert 'length' to numeric. Augur expects numeric types for numeric comparisons.
        if 'length' in metadata.columns:
            log_with_context(logger, "INFO", "Converting 'length' column to numeric.", 
                             {"original_dtype": str(metadata['length'].dtype)})
            metadata['length'] = pd.to_numeric(metadata['length'], errors='coerce')
            # After this, 'length' column will be float type if NaNs were introduced by coercion, 
            # or int if all were valid integers and no NaNs.
            # Records that couldn't be converted will have NaN in 'length'.
            # augur filter should handle NaNs correctly in numeric comparisons (e.g., NaN < X evaluates to False).
            log_with_context(logger, "INFO", "Finished converting 'length' column to numeric.",
                             {"new_dtype": str(metadata['length'].dtype), 
                              "nan_count": int(metadata['length'].isna().sum())}) # Cast nan_count to int for JSON serialization

        # Fill missing values in object/string columns with empty strings.
        # For numeric columns (like 'length' which is now float/int), 
        # NaNs will be handled by to_csv (written as empty strings by default if na_rep is not specified).
        for col in metadata.select_dtypes(include=['object']).columns:
            metadata[col] = metadata[col].fillna("")
        
        # Validate metadata
        validation_result = validate_metadata(metadata, schema)
        
        # Write validation report if path provided
        if args.report:
            os.makedirs(os.path.dirname(args.report), exist_ok=True)
            with open(args.report, 'w') as f:
                json.dump(validation_result, f, indent=2)
        
        # If strict validation is enabled and there are invalid records, exit
        if args.strict and validation_result.get("invalid_records", 0) > 0:
            log_with_context(logger, "ERROR", "Strict validation failed, exiting")
            log_execution_stats(logger, start_time, {
                "operation": "metadata_preparation", 
                "records": len(metadata),
                "invalid_records": validation_result["invalid_records"]
            }, status="validation_failed")
            return 1
        
        # Write standardized metadata
        metadata.to_csv(args.output, sep='\t', index=False)
        
        log_with_context(logger, "INFO", f"Standardized metadata written to {args.output}", {
            "records": len(metadata),
            "columns": list(metadata.columns)
        })
        
        # Log execution statistics
        log_execution_stats(logger, start_time, {
            "operation": "metadata_preparation",
            "records": len(metadata),
            "invalid_records": validation_result.get("invalid_records", 0)
        })
        
        return 0
    
    except Exception as e:
        log_with_context(logger, "CRITICAL", f"Unhandled exception: {e}")
        log_execution_stats(logger, start_time, {"operation": "metadata_preparation", "error": str(e)}, status="failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())