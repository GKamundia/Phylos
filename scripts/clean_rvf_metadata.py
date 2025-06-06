#!/usr/bin/env python3
"""
Clean and enhance RVF metadata for Nextstrain visualization.

This script implements comprehensive data cleaning for RVF sequences:
1. Date standardization to YYYY-MM-DD format
2. Geographic location parsing and country extraction  
3. Latitude/longitude coordinate assignment
4. Data quality validation and reporting
"""

import pandas as pd
import re
import sys
from datetime import datetime
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RVFMetadataCleaner:
    def __init__(self, lat_longs_file):
        """Initialize the metadata cleaner with coordinate data."""
        self.lat_longs = pd.read_csv(lat_longs_file, sep='\t')
        self.country_coords = dict(zip(
            self.lat_longs['location'], 
            zip(self.lat_longs['latitude'], self.lat_longs['longitude'])
        ))
        logger.info(f"Loaded coordinates for {len(self.country_coords)} countries")
        
    def standardize_date(self, date_str):
        """Convert various date formats to YYYY-MM-DD."""
        if pd.isna(date_str) or str(date_str).strip() == '':
            return None
            
        date_str = str(date_str).strip()
        
        # Define date patterns and their corresponding formats
        patterns = [
            # DD-MMM-YYYY (e.g., "11-Nov-2022")
            (r'(\d{1,2})-([A-Za-z]{3})-(\d{4})', '%d-%b-%Y'),
            # DD-MM-YYYY
            (r'(\d{1,2})-(\d{1,2})-(\d{4})', '%d-%m-%Y'),
            # YYYY-MM-DD (already correct)
            (r'(\d{4})-(\d{1,2})-(\d{1,2})', '%Y-%m-%d'),
            # YYYY only
            (r'^(\d{4})$', '%Y'),
            # MM/DD/YYYY
            (r'(\d{1,2})/(\d{1,2})/(\d{4})', '%m/%d/%Y'),
            # YYYY.MM.DD
            (r'(\d{4})\.(\d{1,2})\.(\d{1,2})', '%Y.%m.%d'),
        ]
        
        for pattern, format_str in patterns:
            if re.match(pattern, date_str):
                try:
                    if format_str == '%Y':
                        # For year only, use January 1st
                        return f"{date_str}-01-01"
                    else:
                        parsed_date = datetime.strptime(date_str, format_str)
                        return parsed_date.strftime('%Y-%m-%d')
                except ValueError:
                    continue
                    
        logger.warning(f"Could not parse date: {date_str}")
        return None
    
    def extract_country_from_location(self, geo_location):
        """Extract country from geographic location string."""
        if pd.isna(geo_location) or str(geo_location).strip() == '':
            return None, None
            
        geo_str = str(geo_location).strip()
        
        # Common patterns for geographic locations
        # Pattern 1: "Country: Location" or "Country:Location"
        if ':' in geo_str:
            parts = geo_str.split(':', 1)
            country = parts[0].strip()
            location = parts[1].strip() if len(parts) > 1 else ''
            return country, location
            
        # Pattern 2: "Country, Location" 
        if ',' in geo_str:
            parts = geo_str.split(',', 1)
            country = parts[0].strip()
            location = parts[1].strip() if len(parts) > 1 else ''
            return country, location
            
        # Pattern 3: Just country name
        # Check if it matches known countries
        if geo_str in self.country_coords:
            return geo_str, ''
            
        # Pattern 4: Try to find known country in the string
        for country in self.country_coords.keys():
            if country.lower() in geo_str.lower():
                location = geo_str.replace(country, '').strip().strip(',').strip(':')
                return country, location
                
        # If no pattern matches, treat as unknown location
        return None, geo_str
    
    def assign_coordinates(self, country):
        """Assign latitude and longitude based on country."""
        if pd.isna(country) or str(country).strip() == '':
            return None, None
            
        country_clean = str(country).strip()
        
        # Direct match
        if country_clean in self.country_coords:
            lat, lon = self.country_coords[country_clean]
            return lat, lon
            
        # Try case-insensitive match
        for known_country, coords in self.country_coords.items():
            if known_country.lower() == country_clean.lower():
                return coords
                
        # Try partial match
        for known_country, coords in self.country_coords.items():
            if known_country.lower() in country_clean.lower() or country_clean.lower() in known_country.lower():
                logger.info(f"Partial match: '{country_clean}' -> '{known_country}'")
                return coords
                
        logger.warning(f"No coordinates found for country: {country_clean}")
        return None, None
    
    def clean_metadata(self, metadata_df):
        """Apply all cleaning operations to the metadata DataFrame."""
        logger.info(f"Starting metadata cleaning for {len(metadata_df)} records")
        
        # Create a copy to avoid modifying original
        df = metadata_df.copy()
        
        # 1. Standardize dates
        logger.info("Standardizing collection dates...")
        df['Collection_Date_Original'] = df['Collection_Date'].copy()
        df['Collection_Date'] = df['Collection_Date'].apply(self.standardize_date)
        
        # Count successful date conversions
        valid_dates = df['Collection_Date'].notna().sum()
        logger.info(f"Successfully standardized {valid_dates}/{len(df)} dates")
        
        # 2. Extract countries and locations from Geo_Location
        logger.info("Extracting countries and locations...")
        location_data = df['Geo_Location'].apply(self.extract_country_from_location)
        df['Country_Extracted'] = [item[0] for item in location_data]
        df['Location_Extracted'] = [item[1] for item in location_data]
        
        # Use extracted country if Country field is empty
        df['Country'] = df.apply(lambda row: 
            row['Country_Extracted'] if pd.isna(row['Country']) or str(row['Country']).strip() == '' 
            else row['Country'], axis=1)
        
        # Count successful country extractions
        valid_countries = df['Country'].notna().sum()
        logger.info(f"Successfully extracted {valid_countries}/{len(df)} countries")
        
        # 3. Assign coordinates
        logger.info("Assigning latitude/longitude coordinates...")
        coords = df['Country'].apply(self.assign_coordinates)
        df['latitude'] = [coord[0] if coord[0] is not None else df.iloc[i]['latitude'] for i, coord in enumerate(coords)]
        df['longitude'] = [coord[1] if coord[1] is not None else df.iloc[i]['longitude'] for i, coord in enumerate(coords)]
        
        # Count successful coordinate assignments
        valid_coords = (df['latitude'].notna() & df['longitude'].notna()).sum()
        logger.info(f"Successfully assigned coordinates to {valid_coords}/{len(df)} records")
        
        # 4. Clean strain field (use Accession if strain is empty)
        df['strain'] = df.apply(lambda row: 
            row['Accession'] if pd.isna(row['strain']) or str(row['strain']).strip() == '' 
            else row['strain'], axis=1)
        
        return df
    
    def generate_quality_report(self, original_df, cleaned_df):
        """Generate a quality control report."""
        report = {
            'total_records': len(cleaned_df),
            'dates': {
                'original_valid': original_df['Collection_Date'].notna().sum(),
                'cleaned_valid': cleaned_df['Collection_Date'].notna().sum(),
                'improvement': cleaned_df['Collection_Date'].notna().sum() - original_df['Collection_Date'].notna().sum()
            },
            'countries': {
                'original_valid': original_df['Country'].notna().sum(),
                'cleaned_valid': cleaned_df['Country'].notna().sum(),
                'improvement': cleaned_df['Country'].notna().sum() - original_df['Country'].notna().sum()
            },
            'coordinates': {
                'original_valid': (original_df['latitude'].notna() & original_df['longitude'].notna()).sum(),
                'cleaned_valid': (cleaned_df['latitude'].notna() & cleaned_df['longitude'].notna()).sum(),
                'improvement': (cleaned_df['latitude'].notna() & cleaned_df['longitude'].notna()).sum() - (original_df['latitude'].notna() & original_df['longitude'].notna()).sum()
            }
        }
        
        return report

def main():
    if len(sys.argv) != 4:
        print("Usage: python clean_rvf_metadata.py <input_metadata.tsv> <lat_longs.tsv> <output_metadata.tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    lat_longs_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Load metadata
    logger.info(f"Loading metadata from {input_file}")
    metadata = pd.read_csv(input_file, sep='\t')
    logger.info(f"Loaded {len(metadata)} records")
    
    # Initialize cleaner
    cleaner = RVFMetadataCleaner(lat_longs_file)
    
    # Clean metadata
    cleaned_metadata = cleaner.clean_metadata(metadata)
    
    # Generate quality report
    report = cleaner.generate_quality_report(metadata, cleaned_metadata)
    
    # Print quality report
    print("\n" + "="*60)
    print("RVF METADATA CLEANING REPORT")
    print("="*60)
    print(f"Total records processed: {report['total_records']}")
    print(f"\nDates:")
    print(f"  Original valid: {report['dates']['original_valid']}")
    print(f"  Cleaned valid: {report['dates']['cleaned_valid']}")
    print(f"  Improvement: +{report['dates']['improvement']}")
    print(f"\nCountries:")
    print(f"  Original valid: {report['countries']['original_valid']}")
    print(f"  Cleaned valid: {report['countries']['cleaned_valid']}")
    print(f"  Improvement: +{report['countries']['improvement']}")
    print(f"\nCoordinates:")
    print(f"  Original valid: {report['coordinates']['original_valid']}")
    print(f"  Cleaned valid: {report['coordinates']['cleaned_valid']}")
    print(f"  Improvement: +{report['coordinates']['improvement']}")
    print("="*60)
    
    # Save cleaned metadata
    logger.info(f"Saving cleaned metadata to {output_file}")
    cleaned_metadata.to_csv(output_file, sep='\t', index=False)
    
    # Save quality report
    report_file = output_file.replace('.tsv', '_quality_report.txt')
    with open(report_file, 'w') as f:
        f.write("RVF METADATA CLEANING REPORT\n")
        f.write("="*60 + "\n")
        f.write(f"Total records processed: {report['total_records']}\n")
        f.write(f"\nDates:\n")
        f.write(f"  Original valid: {report['dates']['original_valid']}\n")
        f.write(f"  Cleaned valid: {report['dates']['cleaned_valid']}\n")
        f.write(f"  Improvement: +{report['dates']['improvement']}\n")
        f.write(f"\nCountries:\n")
        f.write(f"  Original valid: {report['countries']['original_valid']}\n")
        f.write(f"  Cleaned valid: {report['countries']['cleaned_valid']}\n")
        f.write(f"  Improvement: +{report['countries']['improvement']}\n")
        f.write(f"\nCoordinates:\n")
        f.write(f"  Original valid: {report['coordinates']['original_valid']}\n")
        f.write(f"  Cleaned valid: {report['coordinates']['cleaned_valid']}\n")
        f.write(f"  Improvement: +{report['coordinates']['improvement']}\n")
    
    logger.info(f"Quality report saved to {report_file}")
    logger.info("Metadata cleaning completed successfully!")

if __name__ == "__main__":
    main()
