#!/usr/bin/env python3
"""
Create reference files for Nextclade from filtered sequences
"""

import os
import pandas as pd
from Bio import SeqIO

# Create reference directory
os.makedirs("nextclade/datasets/rvf", exist_ok=True)

# Read metadata to find the best sequence for each segment
metadata_path = "results/filtered/rvf_metadata.tsv"
sequences_path = "results/filtered/rvf_filtered.fasta"

print(f"Reading metadata from {metadata_path}")
metadata = pd.read_csv(metadata_path, sep='\t')

# Get sequences by segment
L_metadata = metadata[metadata['Segment'] == 'L'].sort_values(by=['Nuc_Completeness', 'Length'], ascending=[False, False])
M_metadata = metadata[metadata['Segment'] == 'M'].sort_values(by=['Nuc_Completeness', 'Length'], ascending=[False, False])
S_metadata = metadata[metadata['Segment'] == 'S'].sort_values(by=['Nuc_Completeness', 'Length'], ascending=[False, False])

# Select reference strains
references = []
if not L_metadata.empty:
    references.append(L_metadata.iloc[0]['strain'])
if not M_metadata.empty:
    references.append(M_metadata.iloc[0]['strain'])
if not S_metadata.empty:
    references.append(S_metadata.iloc[0]['strain'])

print(f"Selected reference strains: {references}")

# Create a reference.fasta file for Nextclade
with open("nextclade/datasets/rvf/reference.fasta", "w") as out_f:
    for record in SeqIO.parse(sequences_path, "fasta"):
        if record.id in references:
            # Get segment from metadata
            segment = metadata[metadata['strain'] == record.id]['Segment'].values[0]
            record.description = f"Rift Valley fever virus {segment} segment reference"
            SeqIO.write(record, out_f, "fasta")
            print(f"Added {record.id} ({segment} segment) to reference file")

# Create a basic QC configuration file
qc_config = """
{
  "schemaVersion": "1.0.0",
  "qc": {
    "missingData": {
      "missingDataThreshold": 1000,
      "scoreBias": 3
    },
    "mixedSites": {
      "mixedSitesThreshold": 10,
      "scoreBias": 3
    },
    "snpClusters": {
      "windowSize": 100,
      "clusterCutOff": 6,
      "scoreWeight": 50
    }
  }
}
"""

with open("nextclade/datasets/rvf/qc.json", "w") as qc_file:
    qc_file.write(qc_config)

print("Created Nextclade reference files in nextclade/datasets/rvf/")
