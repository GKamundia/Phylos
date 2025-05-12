import os
import subprocess
import sys

import snakemake

# Get parameters from Snakemake
output_sequences = snakemake.output.sequences
output_metadata = snakemake.output.metadata
email = snakemake.params.email
search_term = snakemake.params.search_term
max_sequences = snakemake.params.max_sequences
segment = snakemake.params.segment
tracking_file = snakemake.params.tracking_file
archive = snakemake.params.archive
pathogen = snakemake.params.pathogen
incremental = snakemake.params.incremental
log_file = snakemake.log[0]

# Build command
cmd = [
    "python", "scripts/download_sequences.py",
    "--search-term", search_term,
    "--max-sequences", str(max_sequences),
    "--output-sequences", output_sequences,
    "--output-metadata", output_metadata,
    "--email", email,
    "--tracking-file", tracking_file,
    "--archive-dir", "data",
    "--pathogen", pathogen
]

# Add conditional arguments
if segment:
    cmd.extend(["--segment", segment])
if archive:
    cmd.append("--archive")
if incremental:
    cmd.append("--incremental")

# Execute command
with open(log_file, "w") as log:
    result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    
sys.exit(result.returncode)