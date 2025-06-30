import pandas as pd

metadata_path = "results/segments/L/filtered/rvf_L_metadata.tsv"

# Read as raw lines to check for misalignment
with open(metadata_path, encoding="utf-8") as f:
    lines = f.readlines()

header_cols = lines[0].strip().split("\t")
expected_cols = len(header_cols)
bad_rows = []
for i, line in enumerate(lines[1:], 2):
    ncols = len(line.strip("\n").split("\t"))
    if ncols != expected_cols:
        bad_rows.append((i, ncols, line.strip()))

if bad_rows:
    print("Rows with wrong number of columns:")
    for row in bad_rows:
        print(f"Line {row[0]}: {row[1]} columns")
else:
    print("All rows have the correct number of columns.")

# Now load with pandas and check for missing/invalid values
df = pd.read_csv(metadata_path, sep="\t", dtype=str)

required_cols = ["strain", "date", "country", "accession"]
for col in required_cols:
    if col not in df.columns:
        print(f"Missing required column: {col}")
    else:
        missing = df[col].isnull() | (df[col].str.strip() == "")
        if missing.any():
            print(f"Missing values in column '{col}':")
            print(df[missing][[col]])

# Check for invalid dates
if "date" in df.columns:
    invalid_dates = pd.to_datetime(df["date"], errors="coerce").isnull() & (df["date"].str.strip() != "")
    if invalid_dates.any():
        print("Invalid date values:")
        print(df[invalid_dates][["date"]])
    else:
        print("All date values are valid ISO dates.")

print("Metadata integrity check complete.")