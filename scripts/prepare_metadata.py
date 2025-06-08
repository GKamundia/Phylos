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
      # Try various formats (enhanced for NCBI data)
    formats = [
        '%Y-%m-%d', '%Y/%m/%d', '%d-%m-%Y', '%d/%m/%Y', 
        '%m-%d-%Y', '%m/%d/%Y', '%d-%b-%Y', '%d %b %Y', 
        '%b %d %Y', '%B %d %Y', '%d %B %Y', '%d-%B-%Y',
        '%Y-%b-%d', '%Y %b %d', '%b-%Y', '%B-%Y'  # Added month-year formats
    ]
    
    for fmt in formats:
        try:
            date_obj = datetime.strptime(date_str, fmt)
            # For month-year formats, return as YYYY-MM
            if fmt in ['%b-%Y', '%B-%Y']:
                return date_obj.strftime('%Y-%m')
            else:
                return date_obj.strftime('%Y-%m-%d')
        except ValueError:
            continue
    
    # Try to extract year from string if formats fail
    year_match = re.search(r'\b(19|20)\d{2}\b', date_str)
    if year_match:
        log_with_context(logger, "WARNING", f"Could only extract year from date: {date_str} -> {year_match.group()}")
        return year_match.group()
    
    # If all formats fail, log and return original
    log_with_context(logger, "WARNING", f"Could not standardize date: {date_str}")
    return date_str

def extract_country_from_geo_location(geo_location):
    """Extract country from Geo_Location field"""
    if not geo_location or pd.isna(geo_location) or str(geo_location).strip() == "":
        return ""
    
    geo_str = str(geo_location).strip()
    
    # If it contains a colon, take the part before the colon as country
    if ":" in geo_str:
        country = geo_str.split(":")[0].strip()
        return country
    
    # Otherwise, assume the whole string is the country
    return geo_str

def extract_country_from_strain_or_title(isolate, title):
    """Extract country information from strain name or title"""
    
    # Define common country patterns found in RVF data
    country_patterns = {
        r'\b(Kenya|Kenyan)\b': 'Kenya',
        r'\b(South Africa|SA)\b': 'South Africa', 
        r'\b(Egypt|Egyptian)\b': 'Egypt',
        r'\b(Saudi Arabia|Saudi)\b': 'Saudi Arabia',
        r'\b(Madagascar)\b': 'Madagascar',
        r'\b(Senegal)\b': 'Senegal',
        r'\b(Mauritania)\b': 'Mauritania',
        r'\b(Tanzania)\b': 'Tanzania',
        r'\b(Sudan)\b': 'Sudan',
        r'\b(Niger)\b': 'Niger',
        r'\b(Botswana)\b': 'Botswana',
        r'\b(Namibia)\b': 'Namibia',
        r'\b(Uganda)\b': 'Uganda',
        r'\b(Chad)\b': 'Chad',
        r'\b(Mali)\b': 'Mali',
        r'\b(Morocco)\b': 'Morocco',
        r'\b(Central African Republic|CAR)\b': 'Central African Republic',
        r'\b(Rwanda)\b': 'Rwanda',
        r'\b(Burundi)\b': 'Burundi',
        r'\b(Ethiopia)\b': 'Ethiopia',
        r'\b(Somalia)\b': 'Somalia',
        r'\b(Djibouti)\b': 'Djibouti',
        r'\b(Eritrea)\b': 'Eritrea',
        r'\b(Gambia)\b': 'Gambia',
        r'\b(Guinea)\b': 'Guinea',
        r'\b(Ivory Coast|Cote d\'Ivoire)\b': 'Ivory Coast',
        r'\b(Liberia)\b': 'Liberia',
        r'\b(Sierra Leone)\b': 'Sierra Leone',
        r'\b(Burkina Faso)\b': 'Burkina Faso',
        r'\b(Ghana)\b': 'Ghana',
        r'\b(Togo)\b': 'Togo',
        r'\b(Benin)\b': 'Benin',
        r'\b(Nigeria)\b': 'Nigeria',
        r'\b(Cameroon)\b': 'Cameroon',
        r'\b(Gabon)\b': 'Gabon',
        r'\b(Republic of the Congo|Congo)\b': 'Republic of the Congo',
        r'\b(Democratic Republic of the Congo|DRC)\b': 'Democratic Republic of the Congo',
        r'\b(Angola)\b': 'Angola',
        r'\b(Zambia)\b': 'Zambia',
        r'\b(Zimbabwe)\b': 'Zimbabwe',
        r'\b(Mozambique)\b': 'Mozambique',
        r'\b(Malawi)\b': 'Malawi',
        r'\b(Lesotho)\b': 'Lesotho',
        r'\b(Swaziland|Eswatini)\b': 'Eswatini'
    }
    
    # Check isolate/strain field first
    if isolate and not pd.isna(isolate) and str(isolate).strip():
        isolate_str = str(isolate).strip()
        for pattern, country in country_patterns.items():
            if re.search(pattern, isolate_str, re.IGNORECASE):
                return country
    
    # Check title field if no match in isolate
    if title and not pd.isna(title) and str(title).strip():
        title_str = str(title).strip()
        for pattern, country in country_patterns.items():
            if re.search(pattern, title_str, re.IGNORECASE):
                return country
    
    return ""

def create_strain_field(row):
    """Create strain field from available data"""
    # Priority order: Isolate, then Accession
    if 'Isolate' in row and row['Isolate'] and not pd.isna(row['Isolate']) and str(row['Isolate']).strip():
        return str(row['Isolate']).strip()
    elif 'Accession' in row and row['Accession'] and not pd.isna(row['Accession']):
        return str(row['Accession']).strip()
    else:
        return ""

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
        
        # Create standard Nextstrain fields from existing data
        log_with_context(logger, "INFO", "Creating standard Nextstrain fields...")
          # Create 'date' field from Collection_Date only (no fallback to Release_Date)
        if 'date' not in metadata.columns:
            if 'Collection_Date' in metadata.columns:
                log_with_context(logger, "INFO", "Creating standardized 'date' field from 'Collection_Date' only")
                metadata['date'] = metadata['Collection_Date'].apply(standardize_date)
            else:
                log_with_context(logger, "WARNING", "No Collection_Date field found, creating empty 'date' column")
                metadata['date'] = ""
        else:
            # Standardize existing date field
            metadata['date'] = metadata['date'].apply(standardize_date)
        
        # Log date standardization results
        if 'date' in metadata.columns:
            date_stats = metadata['date'].str.len().value_counts()
            log_with_context(logger, "INFO", "Date standardization complete", {
                "complete_dates": int(date_stats.get(10, 0)),  # YYYY-MM-DD (length 10)
                "partial_dates": {
                    "year_month": int(date_stats.get(7, 0)),  # YYYY-MM (length 7)
                    "year": int(date_stats.get(4, 0))  # YYYY (length 4)
                },
                "empty_dates": int(metadata['date'].isna().sum())
            })
        
        # Create 'country' field if missing or enhance existing one
        if 'country' not in metadata.columns:
            if 'Country' in metadata.columns and not metadata['Country'].isna().all():
                metadata['country'] = metadata['Country']
            elif 'Geo_Location' in metadata.columns:
                log_with_context(logger, "INFO", "Extracting country from Geo_Location field")
                metadata['country'] = metadata['Geo_Location'].apply(extract_country_from_geo_location)
                # Also populate the Country field for consistency
                metadata['Country'] = metadata['country']
            else:
                log_with_context(logger, "WARNING", "No geographic information found in standard fields")
                metadata['country'] = ""
                metadata['Country'] = ""
        
        # Enhance country extraction from strain names and titles if countries are still missing
        if 'country' in metadata.columns:
            missing_countries = metadata['country'].isna() | (metadata['country'] == "")
            countries_to_extract = metadata[missing_countries]
            
            if len(countries_to_extract) > 0:
                log_with_context(logger, "INFO", f"Attempting to extract countries from strain names and titles for {len(countries_to_extract)} records")
                
                # Extract countries from isolate and title fields
                for idx, row in countries_to_extract.iterrows():
                    isolate = row.get('Isolate', '')
                    title = row.get('GenBank_Title', '')
                    
                    extracted_country = extract_country_from_strain_or_title(isolate, title)
                    if extracted_country:
                        metadata.at[idx, 'country'] = extracted_country
                        if 'Country' in metadata.columns:
                            metadata.at[idx, 'Country'] = extracted_country
                
                # Log results
                final_missing = metadata['country'].isna() | (metadata['country'] == "")
                extracted_count = len(countries_to_extract) - final_missing.sum()
                log_with_context(logger, "INFO", f"Successfully extracted countries for {extracted_count} additional records")
        
        # Create 'strain' field if missing
        if 'strain' not in metadata.columns:
            log_with_context(logger, "INFO", "Creating strain field from available data")
            metadata['strain'] = metadata.apply(create_strain_field, axis=1)
          # Add virus field if missing (standard for RVF)
        if 'virus' not in metadata.columns:
            metadata['virus'] = 'rvf'
        
        # Create lowercase 'accession' field from 'Accession' if needed for schema compliance
        if 'accession' not in metadata.columns and 'Accession' in metadata.columns:
            metadata['accession'] = metadata['Accession']
        
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