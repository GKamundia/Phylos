import csv

# Input files
minimal_file = "results/segments/L/filtered/rvf_L_metadata_minimal.tsv"
full_file = "results/segments/L/filtered/rvf_L_metadata.tsv"
output_file = "results/segments/L/filtered/rvf_L_metadata_minimal_latlong.tsv"

# Read latitude/longitude from full metadata
latlong_map = {}
with open(full_file, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        strain = row.get('strain') or row.get('accession') or row.get('Accession')
        lat = row.get('latitude', '')
        lon = row.get('longitude', '')
        if strain:
            latlong_map[strain] = (lat, lon)

# Read minimal metadata and add lat/long
with open(minimal_file, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)
    fieldnames = reader.fieldnames + ['latitude', 'longitude']

# Write output
with open(output_file, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for row in rows:
        strain = row['strain']
        lat, lon = latlong_map.get(strain, ('', ''))
        row['latitude'] = lat
        row['longitude'] = lon
        writer.writerow(row)

print(f"Wrote: {output_file}")
