"""
Snakefile for Pathogen-Agnostic Nextstrain pipeline
This Snakefile supports multiple pathogens through the master_config.yaml
"""

import os
import sys
import yaml
from pathlib import Path
from snakemake.utils import min_version

# Ensure minimum Snakemake version for compatibility
min_version("6.0.0")

# Function for safely loading YAML with error handling
def load_yaml(yaml_file):
    try:
        with open(yaml_file) as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading {yaml_file}: {e}", file=sys.stderr)
        sys.exit(1)

# Load master configuration
master_config = load_yaml("config/master_config.yaml")

# Determine active pathogen
active_pathogen = master_config["active_pathogen"]
try:
    pathogen_config = master_config["pathogens"][active_pathogen]
except KeyError:
    print(f"Error: Pathogen '{active_pathogen}' not found in master_config.yaml", file=sys.stderr)
    sys.exit(1)

# Ensure config path exists
if not os.path.exists(pathogen_config["config_path"]):
    print(f"Error: Pathogen config file {pathogen_config['config_path']} not found", file=sys.stderr)
    sys.exit(1)

# Load pathogen-specific configuration
pathogen_specific_config = load_yaml(pathogen_config["config_path"])

# Merge configurations for use in the pipeline
config = {
    "pathogen": active_pathogen,
    "pathogen_name": pathogen_config["name"],
    "has_segments": pathogen_config["has_segments"],
    "data": pathogen_specific_config.get("data", {}),
    "filter": {**master_config["common"]["filter"], **pathogen_specific_config.get("filter", {})},
    "subsample": {**master_config["common"]["subsample"], **pathogen_specific_config.get("subsample", {})},
    "update": master_config["common"]["update"],
    "qc": master_config["common"].get("qc", {}),
    "resources": pathogen_specific_config.get("resources", master_config["common"].get("resources", {}))
}

# Define output prefix based on active pathogen
output_prefix = active_pathogen

# Handle segment configuration
segments = ["L", "M", "S"] if config["has_segments"] and config["data"].get("segment") == "all" else [config["data"].get("segment", "")]
segment_mode = "multi" if len(segments) > 1 and segments[0] != "" else "single"

# Create directories if they don't exist
required_dirs = [
    "logs", 
    "results",
    Path("data/sequences/raw"),
    Path("data/metadata/raw"),
    "results/nextclade",
    "results/qc_reports"
]

for directory in required_dirs:
    os.makedirs(directory, exist_ok=True)

# Create segment-specific directories if needed
if segment_mode == "multi":
    for segment in segments:
        os.makedirs(f"results/segments/{segment}", exist_ok=True)

# Define final output target
rule all:
    input:
        # Standard output for single segment mode
        auspice_json = f"results/auspice/{output_prefix}.json" if segment_mode == "single" else [],
        qc_report = f"results/qc_reports/{output_prefix}_qc_report.html" if segment_mode == "single" else [],
        # Multi-segment outputs if applicable
        segment_jsons = [f"results/segments/{segment}/auspice/{output_prefix}_{segment}.json" for segment in segments] if segment_mode == "multi" else [],
        segment_qc = [f"results/segments/{segment}/qc_reports/{output_prefix}_{segment}_qc_report.html" for segment in segments] if segment_mode == "multi" else [],
        # Combined multi-segment visualization (if applicable)
        combined_json = f"results/auspice/{output_prefix}_combined.json" if segment_mode == "multi" else []

# Include appropriate workflow based on segment mode
include: "workflow/common_rules.smk"

if segment_mode == "multi":
    include: "workflow/multi_segment_rules.smk"

# Download data - segment-aware
rule download_data:
    output:
        sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta",
        metadata = f"data/metadata/raw/{output_prefix}_metadata.tsv"
    params:
        email = config.get("email", "your.email@example.com"),
        search_term = config["data"]["search_term"],
        max_sequences = config["data"]["max_sequences"],
        segment = "all" if segment_mode == "multi" else config["data"].get("segment", ""),
        tracking_file = "config/data_tracking.json",
        archive = config["update"].get("archive", {}).get("enabled", False),
        pathogen = config["pathogen"],
        incremental = config["update"].get("incremental", False)
    log:
        f"logs/download_data_{output_prefix}.log"
    benchmark:
        f"benchmarks/download_data_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("download", {}).get("mem_mb", 2000),
        runtime = config["resources"].get("download", {}).get("runtime", 60)
    script:
        "scripts/run_download.py"

# Clean, validate and standardize metadata
rule prepare_metadata:
    input:
        metadata = f"data/metadata/raw/{output_prefix}_metadata.tsv",
        schema = "config/metadata_schema.json",
        lat_longs = config.get("lat_longs_path", "config/lat_longs.tsv")
    output:
        metadata = f"data/metadata/{output_prefix}_metadata.tsv",
        report = f"data/metadata/{output_prefix}_validation_report.json"
    params:
        strict = False  # Could be parameterized in config
    log:
        f"logs/prepare_metadata_{output_prefix}.log"
    benchmark:
        f"benchmarks/prepare_metadata_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("metadata", {}).get("mem_mb", 2000),
        runtime = config["resources"].get("metadata", {}).get("runtime", 30)
    shell:
        """
        python scripts/prepare_metadata.py \
            {input.metadata} \
            {output.metadata} \
            --schema {input.schema} \
            --lat-longs {input.lat_longs} \
            --report {output.report} \
            {params.strict and '--strict' or ''} \
            > {log} 2>&1
        """

# Filter sequences
rule filter:
    input:
        sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta",
        metadata = f"data/metadata/{output_prefix}_metadata.tsv"
    output:
        sequences = f"results/filtered/{output_prefix}_filtered.fasta",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv"
    params:
        min_length = config["filter"].get("min_length", 0),
        max_n = config["filter"].get("max_n", 100),
        exclude_ids = config["filter"].get("exclude_ids", []),
        exclude_where = lambda w: " ".join([f"--exclude-where \"{condition}\"" for condition in config["filter"].get("exclude_where", [])])
    log:
        f"logs/filter_{output_prefix}.log"
    benchmark:
        f"benchmarks/filter_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("filter", {}).get("mem_mb", 2000)
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --output {output.sequences} \
            --exclude-where "length < {params.min_length}" \
            --exclude-ambiguous-nucleotides-threshold {params.max_n} \
            {params.exclude_where} \
            {'--exclude-ids ' + ' '.join(params.exclude_ids) if params.exclude_ids else ''} \
            --output-metadata {output.metadata} \
            > {log} 2>&1
        """

# Align sequences
rule align:
    input:
        sequences = f"results/filtered/{output_prefix}_filtered.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/filtered/{output_prefix}_{wildcards.segment}_filtered.fasta"
    output:
        alignment = f"results/aligned/{output_prefix}_aligned.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/aligned/{output_prefix}_{wildcards.segment}_aligned.fasta"
    params:
        reference = lambda wildcards: config["data"].get("reference_sequence", {}).get(wildcards.segment if segment_mode == "multi" else config["data"].get("segment", ""), "")
    threads: config["resources"].get("align", {}).get("threads", 4)
    log:
        f"logs/align_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/align_{output_prefix}_{wildcards.segment}.log"
    benchmark:
        f"benchmarks/align_{output_prefix}.txt" if segment_mode == "single" else lambda wildcards: f"benchmarks/align_{output_prefix}_{wildcards.segment}.txt"
    resources:
        mem_mb = config["resources"].get("align", {}).get("mem_mb", 4000)
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --output {output.alignment} \
            {f'--reference-sequence {params.reference}' if params.reference else ''} \
            --nthreads {threads} \
            > {log} 2>&1
        """

# Build tree
rule tree:
    input:
        alignment = f"results/aligned/{output_prefix}_aligned.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/aligned/{output_prefix}_{wildcards.segment}_aligned.fasta"
    output:
        tree = f"results/tree/{output_prefix}_tree.nwk" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/tree/{output_prefix}_{wildcards.segment}_tree.nwk"
    threads: config["resources"].get("tree", {}).get("threads", 4)
    log:
        f"logs/tree_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/tree_{output_prefix}_{wildcards.segment}.log"
    benchmark:
        f"benchmarks/tree_{output_prefix}.txt" if segment_mode == "single" else lambda wildcards: f"benchmarks/tree_{output_prefix}_{wildcards.segment}.txt"
    resources:
        mem_mb = config["resources"].get("tree", {}).get("mem_mb", 8000)
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads} \
            > {log} 2>&1
        """

# Refine tree
rule refine:
    input:
        tree = f"results/tree/{output_prefix}_tree.nwk" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/tree/{output_prefix}_{wildcards.segment}_tree.nwk",
        alignment = f"results/aligned/{output_prefix}_aligned.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/aligned/{output_prefix}_{wildcards.segment}_aligned.fasta",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/filtered/{output_prefix}_{wildcards.segment}_metadata.tsv"
    output:
        tree = f"results/tree/{output_prefix}_refined.nwk" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/tree/{output_prefix}_{wildcards.segment}_refined.nwk",
        node_data = f"results/node_data/{output_prefix}_branch_lengths.json" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/node_data/{output_prefix}_{wildcards.segment}_branch_lengths.json"
    params:
        coalescent = config["refine"].get("coalescent", "opt"),
        clock_rate = config["refine"].get("clock_rate", "")
    threads: config["resources"].get("refine", {}).get("threads", 2)
    log:
        f"logs/refine_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/refine_{output_prefix}_{wildcards.segment}.log"
    benchmark:
        f"benchmarks/refine_{output_prefix}.txt" if segment_mode == "single" else lambda wildcards: f"benchmarks/refine_{output_prefix}_{wildcards.segment}.txt"
    resources:
        mem_mb = config["resources"].get("refine", {}).get("mem_mb", 4000)
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            {f'--clock-rate {params.clock_rate}' if params.clock_rate else ''} \
            --date-confidence \
            > {log} 2>&1
        """