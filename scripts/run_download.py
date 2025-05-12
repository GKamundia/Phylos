import os
import sys
import json
import time
from datetime import datetime
from Bio import Entrez, SeqIO
#import snakemake

# Ensure an email is provided for Entrez queries
entrez_email = snakemake.params.get("email", "gkamush69@gmai.com")
if entrez_email == "gkamush69@gmai.com" or not entrez_email:
    print(
        "Warning: Using default NCBI Entrez email or no email provided. "
        "Please set a proper email in your config (e.g., config.yaml under 'email' or within params for the rule).",
        file=sys.stderr
    )
Entrez.email = entrez_email

# Access output file paths
output_sequences_path = snakemake.output.sequences
output_metadata_path = snakemake.output.metadata

# Access log file path
log_file_path = snakemake.log[0]

# Access params
search_term = snakemake.params.search_term
max_sequences = int(snakemake.params.max_sequences)
segment = snakemake.params.segment
tracking_file = snakemake.params.tracking_file
archive = snakemake.params.archive
pathogen = snakemake.params.pathogen
incremental = snakemake.params.incremental

print(f"Starting data download: {datetime.now()}")
print(f"Output sequences file: {output_sequences_path}")
print(f"Output metadata file: {output_metadata_path}")
print(f"Search term: {search_term}")
print(f"Max sequences: {max_sequences}")

# Placeholder for your actual download logic
# Ensure you use the variables like output_sequences_path, output_metadata_path,
# search_term, max_sequences, etc., derived from the snakemake object.

# Example: Creating dummy output files for testing
# Remove this section when you have your actual download logic
with open(output_sequences_path, "w") as f_seq:
    f_seq.write(">dummy_seq\nACGT\n")
with open(output_metadata_path, "w") as f_meta:
    f_meta.write("strain\tdate\n")
    f_meta.write("dummy_strain\t2023-01-01\n")

print(f"Data download script finished: {datetime.now()}")