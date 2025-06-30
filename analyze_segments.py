import pandas as pd

df = pd.read_csv("results/segments/L/filtered/rvf_L_metadata.tsv", sep="\t")
for col in ["Length", "latitude", "longitude"]:
    non_numeric = df[~df[col].astype(str).str.replace('.', '', 1).str.replace('-', '', 1).str.isnumeric()]
    if not non_numeric.empty:
        print(f"Non-numeric values in {col}:")
        print(non_numeric[[col]])