"""
Script to generate rvf_M_metadata_latlong.tsv from rvf_M_metadata.tsv
Extracts columns: strain, country, latitude, longitude for all entries.
Sets the strain column to be identical to the Accession column (GenBank accession).
"""
import pandas as pd

input_file = "results/segments/M/filtered/rvf_M_metadata.tsv"
output_file = "results/segments/M/filtered/rvf_M_metadata_latlong.tsv"

df = pd.read_csv(input_file, sep='\t', dtype=str)
# Set the strain column to be identical to the Accession column
df['strain'] = df['Accession']
# Only keep rows where all required fields are present and non-empty
required = ['strain', 'country', 'latitude', 'longitude']
df_out = df.dropna(subset=required)
df_out = df_out[["strain", "country", "latitude", "longitude"]]
df_out.to_csv(output_file, sep='\t', index=False)
print(f"Saved file with columns: strain, country, latitude, longitude to {output_file}")
