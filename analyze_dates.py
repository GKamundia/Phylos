#!/usr/bin/env python3
import pandas as pd

# Load the prepared metadata
df = pd.read_csv('results/test_metadata.tsv', sep='\t')

print('=== COMPREHENSIVE METADATA ANALYSIS ===')
print(f'Total records: {len(df)}')
print()

print('=== DATE ANALYSIS ===')
print(f'Records with Collection_Date: {df["Collection_Date"].notna().sum()}')
print(f'Records with Release_Date: {df["Release_Date"].notna().sum()}')
print(f'Records with standardized date: {(df["date"] != "").sum()}')
print(f'Empty dates: {(df["date"] == "").sum()}')
print()

print('Date field length distribution:')
print(df['date'].str.len().value_counts().sort_index())
print()

print('Records using Release_Date fallback (no Collection_Date):')
empty_collection = df[df['Collection_Date'].isna() | (df['Collection_Date'] == '')]
print(f'Count: {len(empty_collection)}')
if len(empty_collection) > 0:
    print(empty_collection[['Accession', 'Collection_Date', 'Release_Date', 'date']].head())
print()

print('=== GEOGRAPHIC COORDINATE ANALYSIS ===')
has_coords = (df['latitude'].notna()) & (df['longitude'].notna())
print(f'Records with coordinates: {has_coords.sum()}')
print(f'Missing coordinates: {(~has_coords).sum()}')
print()

print('Countries with coordinates:')
countries_with_coords = df[has_coords]['country'].value_counts()
print(countries_with_coords.head(10))
print()

print('Countries missing coordinates:')
missing_coords = df[~has_coords]['country'].value_counts()
if len(missing_coords) > 0:
    print(missing_coords)
else:
    print('All countries have coordinates!')
print()

print('=== VALIDATION SUMMARY ===')
print(f'Total records processed: {len(df)}')
print(f'Records with coordinates: {has_coords.sum()} ({has_coords.sum()/len(df)*100:.1f}%)')
print(f'Complete dates (YYYY-MM-DD): {(df["date"].str.len() == 10).sum()} ({(df["date"].str.len() == 10).sum()/len(df)*100:.1f}%)')
print(f'Year-only dates: {(df["date"].str.len() == 4).sum()} ({(df["date"].str.len() == 4).sum()/len(df)*100:.1f}%)')
print(f'Records using Collection_Date: {df["Collection_Date"].notna().sum()} ({df["Collection_Date"].notna().sum()/len(df)*100:.1f}%)')
print(f'Records using Release_Date fallback: {len(empty_collection)} ({len(empty_collection)/len(df)*100:.1f}%)')
