#!/usr/bin/env python3
"""
Final validation script for enhanced RVF Nextstrain pipeline
Tests the complete data cleaning integration and verifies results
"""

import pandas as pd
import json
import os
from pathlib import Path

def validate_enhanced_pipeline():
    """Validate the enhanced pipeline integration and results"""
    
    print("🔍 Enhanced RVF Nextstrain Pipeline Validation")
    print("=" * 50)
      # Define paths
    base_dir = Path(".")
    filtered_metadata = base_dir / "results/filtered/rvf_metadata.tsv"
    enhanced_metadata = base_dir / "data/metadata/rvf_metadata.tsv"
    
    validation_results = {
        "pipeline_integration": {},
        "data_quality": {},
        "geographic_enhancement": {},
        "field_standardization": {}
    }
    
    # Test 1: Pipeline Integration
    print("\n📊 Testing Pipeline Integration...")
    
    if enhanced_metadata.exists():
        enhanced_df = pd.read_csv(enhanced_metadata, sep='\t', dtype=str)
        validation_results["pipeline_integration"]["enhanced_metadata_exists"] = True
        validation_results["pipeline_integration"]["enhanced_records"] = len(enhanced_df)
        print(f"   ✅ Enhanced metadata exists: {len(enhanced_df)} records")
    else:
        validation_results["pipeline_integration"]["enhanced_metadata_exists"] = False
        print("   ❌ Enhanced metadata not found")
        return validation_results
    
    if filtered_metadata.exists():
        filtered_df = pd.read_csv(filtered_metadata, sep='\t', dtype=str)
        validation_results["pipeline_integration"]["filtered_metadata_exists"] = True
        validation_results["pipeline_integration"]["filtered_records"] = len(filtered_df)
        print(f"   ✅ Filtered metadata exists: {len(filtered_df)} records")
    else:
        validation_results["pipeline_integration"]["filtered_metadata_exists"] = False
        print("   ❌ Filtered metadata not found")
        return validation_results
    
    # Test 2: Data Quality Improvements
    print("\n🎯 Testing Data Quality Improvements...")
    
    # Check date standardization
    if 'date' in enhanced_df.columns:
        date_formats = enhanced_df['date'].dropna()
        complete_dates = date_formats[date_formats.str.match(r'^\d{4}-\d{2}-\d{2}$', na=False)]
        year_only = date_formats[date_formats.str.match(r'^\d{4}$', na=False)]
        
        validation_results["data_quality"]["total_dates"] = len(date_formats)
        validation_results["data_quality"]["complete_dates"] = len(complete_dates)
        validation_results["data_quality"]["year_only_dates"] = len(year_only)
        
        print(f"   ✅ Date standardization: {len(complete_dates)} complete, {len(year_only)} year-only")
    
    # Check country extraction
    if 'country' in enhanced_df.columns:
        countries_present = enhanced_df['country'].dropna()
        countries_with_values = countries_present[countries_present != '']
        
        validation_results["data_quality"]["total_countries"] = len(countries_with_values)
        validation_results["data_quality"]["unique_countries"] = len(countries_with_values.unique())
        
        print(f"   ✅ Country extraction: {len(countries_with_values)} records, {len(countries_with_values.unique())} unique countries")
        print(f"      Top countries: {list(countries_with_values.value_counts().head(3).index)}")
    
    # Test 3: Geographic Enhancement
    print("\n🌍 Testing Geographic Enhancement...")
    
    # Check coordinates in enhanced metadata
    lat_cols = [col for col in enhanced_df.columns if 'latitude' in col.lower()]
    lon_cols = [col for col in enhanced_df.columns if 'longitude' in col.lower()]
    
    if lat_cols and lon_cols:
        lat_col, lon_col = lat_cols[0], lon_cols[0]
        coords_present = enhanced_df[
            (enhanced_df[lat_col].notna()) & 
            (enhanced_df[lat_col] != '') &
            (enhanced_df[lon_col].notna()) & 
            (enhanced_df[lon_col] != '')
        ]
        
        validation_results["geographic_enhancement"]["enhanced_with_coords"] = len(coords_present)
        print(f"   ✅ Enhanced metadata with coordinates: {len(coords_present)} records")
    
    # Check coordinates in filtered metadata
    lat_cols_filtered = [col for col in filtered_df.columns if 'latitude' in col.lower()]
    lon_cols_filtered = [col for col in filtered_df.columns if 'longitude' in col.lower()]
    
    if lat_cols_filtered and lon_cols_filtered:
        lat_col_f, lon_col_f = lat_cols_filtered[0], lon_cols_filtered[0]
        coords_filtered = filtered_df[
            (filtered_df[lat_col_f].notna()) & 
            (filtered_df[lat_col_f] != '') &
            (filtered_df[lon_col_f].notna()) & 
            (filtered_df[lon_col_f] != '')
        ]
        
        validation_results["geographic_enhancement"]["filtered_with_coords"] = len(coords_filtered)
        print(f"   ✅ Filtered metadata with coordinates: {len(coords_filtered)} records")
        
        if len(coords_filtered) > 0:
            print("   📍 Example coordinates:")
            for idx, row in coords_filtered.head(3).iterrows():
                country = row.get('country', 'Unknown')
                strain = row.get('strain', 'Unknown')
                lat = row[lat_col_f]
                lon = row[lon_col_f] 
                print(f"      {country} | {strain} | {lat}, {lon}")
    
    # Test 4: Field Standardization
    print("\n📋 Testing Field Standardization...")
    
    required_fields = ['strain', 'virus', 'date', 'country']
    field_status = {}
    
    for field in required_fields:
        if field in enhanced_df.columns:
            non_empty = enhanced_df[field].dropna()
            non_empty = non_empty[non_empty != '']
            field_status[field] = len(non_empty)
            print(f"   ✅ {field}: {len(non_empty)} records with data")
        else:
            field_status[field] = 0
            print(f"   ❌ {field}: Missing field")
    
    validation_results["field_standardization"] = field_status
    
    # Test 5: Integration Verification
    print("\n🔗 Testing Integration Verification...")
    
    # Check if enhanced cleaning results propagated to filtered data
    if 'country' in filtered_df.columns:
        filtered_countries = filtered_df['country'].dropna()
        filtered_countries = filtered_countries[filtered_countries != '']
        print(f"   ✅ Countries in filtered data: {len(filtered_countries)} records")
        
        if len(filtered_countries) > 0:
            print(f"   🌍 Countries represented: {list(filtered_countries.unique())}")
    
    # Summary
    print("\n" + "=" * 50)
    print("📊 VALIDATION SUMMARY")
    print("=" * 50)
    
    total_tests = 0
    passed_tests = 0
    
    # Count test results
    if validation_results["pipeline_integration"].get("enhanced_metadata_exists"):
        passed_tests += 1
    total_tests += 1
    
    if validation_results["pipeline_integration"].get("filtered_metadata_exists"):
        passed_tests += 1
    total_tests += 1
    
    if validation_results["data_quality"].get("total_countries", 0) > 0:
        passed_tests += 1
    total_tests += 1
    
    if validation_results["geographic_enhancement"].get("filtered_with_coords", 0) > 0:
        passed_tests += 1
    total_tests += 1
    
    if all(validation_results["field_standardization"].get(field, 0) > 0 for field in required_fields):
        passed_tests += 1
    total_tests += 1
    
    print(f"✅ Tests Passed: {passed_tests}/{total_tests}")
    
    if passed_tests == total_tests:
        print("🎉 ALL TESTS PASSED - Enhanced pipeline integration successful!")
    else:
        print(f"⚠️  {total_tests - passed_tests} tests failed - Review integration")
    
    # Save validation results
    results_file = base_dir / "enhanced_validation_results.json"
    with open(results_file, 'w') as f:
        json.dump(validation_results, f, indent=2)
    print(f"\n📄 Detailed results saved to: {results_file}")
    
    return validation_results

if __name__ == "__main__":
    validate_enhanced_pipeline()
