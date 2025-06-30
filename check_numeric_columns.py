import pandas as pd

# Path to your metadata file
metadata_path = "results/segments/L/filtered/rvf_L_metadata.tsv"

# Columns to check
numeric_columns = ["Length", "latitude", "longitude"]

# Read the TSV file
df = pd.read_csv(metadata_path, sep="\t", dtype=str)

for col in numeric_columns:
    # Remove whitespace and check if value is numeric (allow negative and decimal)
    non_numeric = df[~df[col].str.strip().replace('', 'NaN').str.replace('.', '', 1).str.replace('-', '', 1).str.isnumeric()]
    if not non_numeric.empty:
        print(f"\nNon-numeric values found in column '{col}':")
        print(non_numeric[[col]])
    else:
        print(f"No non-numeric values found in column '{col}'.")

print("\nCheck complete.")