#!/usr/bin/env python3
import csv
from collections import Counter

metadata_path = 'results/segments/L/filtered/rvf_L_metadata.tsv'

date_counts = Counter()

with open(metadata_path, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    try:
        date_idx = header.index('date')
    except ValueError:
        print("No 'date' column found.")
        exit(1)
    total = 0
    valid = 0
    valid_dates = []
    invalid_dates = []
    for row in reader:
        if len(row) > date_idx:
            date = row[date_idx]
            if date and (len(date) == 10 or len(date) == 4):
                valid += 1
                valid_dates.append(date)
                date_counts[date] += 1
            else:
                invalid_dates.append(date)
            total += 1
    print(f"Valid date entries: {valid}")
    print(f"Total sequences: {total}")
    if total > 0:
        pct = round((valid / total) * 100, 1)
        print(f"Percent valid: {pct}%")
    print("\nSample of valid dates:", valid_dates[:5])
    print("Sample of invalid dates:", invalid_dates[:5])
    print("\n=== Unique date distribution (top 10) ===")
    for date, count in date_counts.most_common(10):
        print(f"{date}: {count}")
    print(f"Total unique dates: {len(date_counts)}")

import pandas as pd

metadata_path = r'c:\Users\Anarchy\Documents\Data_Science\NextStrain\rvf-nextstrain\results\segments\L\filtered\rvf_L_metadata.tsv'

metadata = pd.read_csv(metadata_path, sep='\t', dtype=str)
if 'accession' in metadata.columns and 'strain' in metadata.columns:
    metadata['strain'] = metadata['accession']
    metadata.to_csv(metadata_path, sep='\t', index=False)
    print("Updated 'strain' column to match 'accession' for all rows.")
else:
    raise ValueError("'accession' and/or 'strain' columns not found in metadata.")