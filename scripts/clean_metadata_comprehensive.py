#!/usr/bin/env python3
"""
Comprehensive metadata cleaning script for RVF Nextstrain pipeline
Addresses data cleaning issues identified in pipeline testing:
1. Date format standardization (YYYY-MM-DD)
2. Geographic location parsing (Country extraction from various fields)
3. Latitude/longitude coordinate integration
4. Strain field creation
5. Missing data handling
"""

import os
import sys
import re
import json
import argparse
import pandas as pd
from datetime import datetime
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Comprehensive metadata cleaning for RVF pipeline")
    parser.add_argument("input", help="Input metadata TSV file")
    parser.add_argument("output", help="Output cleaned metadata TSV file")
    parser.add_argument("--lat-longs", help="TSV file with latitude/longitude data", 
                       default="config/lat_longs.tsv")
    parser.add_argument("--report", help="Output path for cleaning report JSON")
    
    return parser.parse_args()

def load_lat_longs(lat_longs_file):
    """Load latitude/longitude reference data"""
    try:
        if not lat_longs_file or not os.path.exists(lat_longs_file):
            logger.warning(f"Lat/long file not found: {lat_longs_file}")
            return {}
        
        # Read lat_longs file
        lat_longs_df = pd.read_csv(lat_longs_file, sep='\t')
        logger.info(f"Loaded {len(lat_longs_df)} geographic coordinate entries")
        
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
        
        return lat_longs
    
    except Exception as e:
        logger.error(f"Error loading lat/longs: {e}")
        return {}

def standardize_date(date_str):
    """Convert various date formats to YYYY-MM-DD format"""
    if not date_str or pd.isna(date_str) or str(date_str).strip() == "":
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
    
    # Try various date formats commonly found in NCBI data
    formats = [
        '%Y-%m-%d', '%Y/%m/%d', '%d-%m-%Y', '%d/%m/%Y', 
        '%m-%d-%Y', '%m/%d/%Y', '%d-%b-%Y', '%d %b %Y', 
        '%b %d %Y', '%B %d %Y', '%d %B %Y', '%d-%B-%Y',
        '%Y-%b-%d', '%Y %b %d'
    ]
    
    for fmt in formats:
        try:
            date_obj = datetime.strptime(date_str, fmt)
            return date_obj.strftime('%Y-%m-%d')
        except ValueError:
            continue
    
    # Try to extract year from string if formats fail
    year_match = re.search(r'\b(19|20)\d{2}\b', date_str)
    if year_match:
        logger.warning(f"Could only extract year from date: {date_str} -> {year_match.group()}")
        return year_match.group()
    
    logger.warning(f"Could not standardize date: {date_str}")
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

def extract_location_from_geo_location(geo_location):
    """Extract specific location from Geo_Location field"""
    if not geo_location or pd.isna(geo_location) or str(geo_location).strip() == "":
        return ""
    
    geo_str = str(geo_location).strip()
    
    # If it contains a colon, take the part after the colon as location
    if ":" in geo_str:
        parts = geo_str.split(":", 1)
        if len(parts) > 1:
            return parts[1].strip()
    
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

def add_geographic_coordinates(metadata_df, lat_longs):
    """Add latitude and longitude coordinates to metadata"""
    if not lat_longs:
        logger.warning("No geographic coordinate data available")
        return metadata_df
    
    # Add lat/long columns if not present
    for col in ['latitude', 'longitude']:
        if col not in metadata_df.columns:
            metadata_df[col] = ""
    
    coords_added = 0
    coords_missing = 0
    
    for idx, row in metadata_df.iterrows():
        coords_found = False
        
        # Try to match by country
        if 'Country' in row and row['Country'] and not pd.isna(row['Country']) and str(row['Country']).strip():
            country = str(row['Country']).strip()
            if 'country' in lat_longs and country in lat_longs['country']:
                lat, long = lat_longs['country'][country]
                metadata_df.at[idx, 'latitude'] = lat
                metadata_df.at[idx, 'longitude'] = long
                coords_found = True
                coords_added += 1
        
        if not coords_found:
            coords_missing += 1
    
    logger.info(f"Added coordinates to {coords_added} records, {coords_missing} records missing coordinates")
    return metadata_df

def create_nextstrain_standard_fields(metadata_df):
    """Create standard Nextstrain fields from existing data"""
    
    # Create 'date' field from Collection_Date
    if 'Collection_Date' in metadata_df.columns:
        logger.info("Creating standardized 'date' field from 'Collection_Date'")
        metadata_df['date'] = metadata_df['Collection_Date'].apply(standardize_date)
    elif 'Release_Date' in metadata_df.columns:
        logger.info("Creating standardized 'date' field from 'Release_Date'")
        metadata_df['date'] = metadata_df['Release_Date'].apply(standardize_date)
    else:
        logger.warning("No date field found, creating empty 'date' column")
        metadata_df['date'] = ""
      # Create 'country' field if missing
    if 'country' not in metadata_df.columns:
        if 'Country' in metadata_df.columns and not metadata_df['Country'].isna().all():
            metadata_df['country'] = metadata_df['Country']
        elif 'Geo_Location' in metadata_df.columns:
            logger.info("Extracting country from Geo_Location field")
            metadata_df['country'] = metadata_df['Geo_Location'].apply(extract_country_from_geo_location)
            # Also populate the Country field for consistency
            metadata_df['Country'] = metadata_df['country']
        else:
            logger.warning("No geographic information found in standard fields")
            metadata_df['country'] = ""
            metadata_df['Country'] = ""
    
    # Enhance country extraction from strain names and titles if countries are still missing
    if 'country' in metadata_df.columns:
        missing_countries = metadata_df['country'].isna() | (metadata_df['country'] == "")
        countries_to_extract = metadata_df[missing_countries]
        
        if len(countries_to_extract) > 0:
            logger.info(f"Attempting to extract countries from strain names and titles for {len(countries_to_extract)} records")
            
            # Extract countries from isolate and title fields
            for idx, row in countries_to_extract.iterrows():
                isolate = row.get('Isolate', '')
                title = row.get('GenBank_Title', '')
                
                extracted_country = extract_country_from_strain_or_title(isolate, title)
                if extracted_country:
                    metadata_df.at[idx, 'country'] = extracted_country
                    metadata_df.at[idx, 'Country'] = extracted_country
            
            # Log results
            final_missing = metadata_df['country'].isna() | (metadata_df['country'] == "")
            extracted_count = len(countries_to_extract) - final_missing.sum()
            logger.info(f"Successfully extracted countries for {extracted_count} additional records")
    
    # Create 'strain' field if missing
    if 'strain' not in metadata_df.columns:
        logger.info("Creating 'strain' field from available data")
        metadata_df['strain'] = metadata_df.apply(create_strain_field, axis=1)
    
    # Create 'virus' field
    if 'virus' not in metadata_df.columns:
        metadata_df['virus'] = 'rvf'  # Standard for RVF virus
    
    # Create 'accession' field (lowercase) if missing
    if 'accession' not in metadata_df.columns and 'Accession' in metadata_df.columns:
        metadata_df['accession'] = metadata_df['Accession']
    
    return metadata_df

def validate_and_clean_data(metadata_df):
    """Perform final validation and cleaning"""
    
    # Convert length to numeric
    if 'Length' in metadata_df.columns:
        logger.info("Converting Length field to numeric")
        metadata_df['length'] = pd.to_numeric(metadata_df['Length'], errors='coerce')
    
    # Fill empty string values for object columns
    for col in metadata_df.select_dtypes(include=['object']).columns:
        metadata_df[col] = metadata_df[col].fillna("")
    
    # Remove any completely empty rows
    initial_count = len(metadata_df)
    metadata_df = metadata_df.dropna(how='all')
    final_count = len(metadata_df)
    
    if initial_count != final_count:
        logger.info(f"Removed {initial_count - final_count} completely empty rows")
    
    return metadata_df

def generate_cleaning_report(original_df, cleaned_df, coords_added):
    """Generate a comprehensive cleaning report"""
    
    report = {
        "timestamp": datetime.now().isoformat(),
        "input_records": len(original_df),
        "output_records": len(cleaned_df),
        "columns_added": [],
        "data_transformations": {
            "dates_standardized": 0,
            "countries_extracted": 0,
            "coordinates_added": coords_added,
            "strains_created": 0
        },
        "data_quality": {
            "records_with_dates": 0,
            "records_with_countries": 0,
            "records_with_coordinates": 0,
            "records_with_strains": 0
        }
    }
    
    # Count standardized dates
    if 'date' in cleaned_df.columns:
        report["data_quality"]["records_with_dates"] = int((cleaned_df['date'] != "").sum())
    
    # Count countries
    if 'country' in cleaned_df.columns:
        report["data_quality"]["records_with_countries"] = int((cleaned_df['country'] != "").sum())
    
    # Count coordinates
    if 'latitude' in cleaned_df.columns:
        report["data_quality"]["records_with_coordinates"] = int((cleaned_df['latitude'] != "").sum())
    
    # Count strains
    if 'strain' in cleaned_df.columns:
        report["data_quality"]["records_with_strains"] = int((cleaned_df['strain'] != "").sum())
    
    # Identify new columns
    original_cols = set(original_df.columns)
    cleaned_cols = set(cleaned_df.columns)
    report["columns_added"] = list(cleaned_cols - original_cols)
    
    return report

def main():
    """Main cleaning function"""
    args = parse_args()
    
    try:
        # Load input metadata
        logger.info(f"Loading metadata from {args.input}")
        metadata_df = pd.read_csv(args.input, sep='\t', dtype=str)
        original_df = metadata_df.copy()
        
        logger.info(f"Loaded {len(metadata_df)} records with {len(metadata_df.columns)} columns")
        
        # Load geographic coordinate reference data
        lat_longs = load_lat_longs(args.lat_longs)
        
        # Create standard Nextstrain fields
        logger.info("Creating standard Nextstrain fields...")
        metadata_df = create_nextstrain_standard_fields(metadata_df)
        
        # Add geographic coordinates
        logger.info("Adding geographic coordinates...")
        coords_added = 0
        if lat_longs:
            metadata_df = add_geographic_coordinates(metadata_df, lat_longs)
            coords_added = int((metadata_df['latitude'] != "").sum())
        
        # Final validation and cleaning
        logger.info("Performing final validation and cleaning...")
        metadata_df = validate_and_clean_data(metadata_df)
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        
        # Write cleaned metadata
        metadata_df.to_csv(args.output, sep='\t', index=False)
        logger.info(f"Cleaned metadata written to {args.output}")
        
        # Generate and save cleaning report
        if args.report:
            report = generate_cleaning_report(original_df, metadata_df, coords_added)
            os.makedirs(os.path.dirname(args.report), exist_ok=True)
            with open(args.report, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Cleaning report written to {args.report}")
            
            # Print summary
            print("\n=== METADATA CLEANING SUMMARY ===")
            print(f"Input records: {report['input_records']}")
            print(f"Output records: {report['output_records']}")
            print(f"Records with dates: {report['data_quality']['records_with_dates']}")
            print(f"Records with countries: {report['data_quality']['records_with_countries']}")
            print(f"Records with coordinates: {report['data_quality']['records_with_coordinates']}")
            print(f"Records with strains: {report['data_quality']['records_with_strains']}")
            print(f"New columns added: {', '.join(report['columns_added'])}")
        
        logger.info("Metadata cleaning completed successfully")
        return 0
        
    except Exception as e:
        logger.error(f"Error during metadata cleaning: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
