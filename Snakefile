"""
Snakefile for Pathogen-Agnostic Nextstrain pipeline
This Snakefile supports multiple pathogens through the master_config.yaml
"""

import os
import sys
import yaml
from pathlib import Path

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
}

# Define output prefix based on active pathogen
output_prefix = active_pathogen

# Create directories if they don't exist
for directory in ["logs", "results", Path("data/sequences/raw"), Path("data/metadata/raw")]:
    os.makedirs(directory, exist_ok=True)

# Define final output target
rule all:
    input:
        auspice_json = f"results/auspice/{output_prefix}.json",
        qc_report = f"results/qc_reports/{output_prefix}_qc_report.html"

# Download data from NCBI with archiving and incremental updates
rule download_data:
    output:
        sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta",
        metadata = f"data/metadata/raw/{output_prefix}_metadata.tsv"
    params:
        email = "your.email@example.com",
        search_term = config["data"]["search_term"],
        max_sequences = config["data"]["max_sequences"],
        segment = config["data"].get("segment", ""),
        tracking_file = "config/data_tracking.json",
        archive = config["update"].get("archive", {}).get("enabled", False),
        pathogen = config["pathogen"],
        incremental = config["update"].get("incremental", False)
    log:
        f"logs/download_data_{output_prefix}.log"
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
    log:
        f"logs/prepare_metadata_{output_prefix}.log"
    shell:
        """
        python scripts/prepare_metadata.py \
            {input.metadata} \
            {output.metadata} \
            --schema {input.schema} \
            --lat-longs {input.lat_longs} \
            --report {output.report} \
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
        max_n = config["filter"].get("max_n", 100)
    log:
        f"logs/filter_{output_prefix}.log"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --output {output.sequences} \
            --exclude-where "length < {params.min_length}" \
            --exclude-ambiguous-nucleotides-threshold {params.max_n} \
            --output-metadata {output.metadata} \
            > {log} 2>&1
        """

# Align sequences
rule align:
    input:
        sequences = f"results/filtered/{output_prefix}_filtered.fasta"
    output:
        alignment = f"results/aligned/{output_prefix}_aligned.fasta"
    threads: 4
    log:
        f"logs/align_{output_prefix}.log"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --output {output.alignment} \
            --nthreads {threads} \
            > {log} 2>&1
        """

# Build tree
rule tree:
    input:
        alignment = f"results/aligned/{output_prefix}_aligned.fasta"
    output:
        tree = f"results/tree/{output_prefix}_tree.nwk"
    threads: 4
    log:
        f"logs/tree_{output_prefix}.log"
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
        tree = f"results/tree/{output_prefix}_tree.nwk",
        alignment = f"results/aligned/{output_prefix}_aligned.fasta",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv"
    output:
        tree = f"results/tree/{output_prefix}_refined.nwk",
        node_data = f"results/node_data/{output_prefix}_branch_lengths.json"
    log:
        f"logs/refine_{output_prefix}.log"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent opt \
            --date-confidence \
            > {log} 2>&1
        """

# Reconstruct ancestral traits (e.g., country)
rule traits:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv"
    output:
        node_data = f"results/node_data/{output_prefix}_traits.json"
    log:
        f"logs/traits_{output_prefix}.log"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns country \
            > {log} 2>&1
        """

# Export to auspice
rule export:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv",
        branch_lengths = f"results/node_data/{output_prefix}_branch_lengths.json",
        traits = f"results/node_data/{output_prefix}_traits.json"
    output:
        auspice_json = f"results/auspice/{output_prefix}.json"
    params:
        auspice_config = pathogen_config.get("auspice_config_path", "config/auspice_config.json"),
        lat_longs = pathogen_config.get("lat_longs_path", "config/lat_longs.tsv")
    log:
        f"logs/export_{output_prefix}.log"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} \
            --auspice-config {params.auspice_config} \
            --lat-longs {params.lat_longs} \
            --output {output.auspice_json} \
            > {log} 2>&1
        """

# Run Nextclade for sequence QC, clade assignment, and outlier detection
rule nextclade_qc:
    input:
        sequences = f"results/filtered/{output_prefix}_filtered.fasta"
    output:
        json = f"results/nextclade/{output_prefix}_nextclade.json",
        tsv = f"results/nextclade/{output_prefix}_nextclade.tsv",
        aligned = f"results/nextclade/{output_prefix}_aligned.fasta",
        passed = f"results/nextclade/{output_prefix}_passed.fasta"
    params:
        segment = config["data"].get("segment", ""),
        dataset_dir = lambda w: f"nextclade/datasets/{config['pathogen']}/",
        min_qc_score = config.get("qc", {}).get("nextclade", {}).get("min_qc_score", 80),
        outdir = "results/nextclade"
    log:
        f"logs/nextclade_qc_{output_prefix}.log"
    benchmark:
        f"benchmarks/nextclade_qc_{output_prefix}.txt"
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p {params.outdir}
        
        # Run Nextclade
        nextclade run \
            --input-dataset {params.dataset_dir} \
            --output-json {output.json} \
            --output-tsv {output.tsv} \
            --output-aligned {output.aligned} \
            --input-fasta {input.sequences} \
            --include-reference \
            > {log} 2>&1
        
        # Filter sequences by QC score and create passed.fasta
        python scripts/filter_nextclade_results.py \
            --input-json {output.json} \
            --input-fasta {input.sequences} \
            --output-fasta {output.passed} \
            --min-qc-score {params.min_qc_score} \
            --segment {params.segment} \
            >> {log} 2>&1
        """

# Generate QC report (after Nextclade and other filters)
rule generate_qc_report:
    input:
        raw_sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta",
        filtered_sequences = f"results/filtered/{output_prefix}_filtered.fasta",
        nextclade_json = f"results/nextclade/{output_prefix}_nextclade.json",
        nextclade_passed = f"results/nextclade/{output_prefix}_passed.fasta",
        raw_metadata = f"data/metadata/raw/{output_prefix}_metadata.tsv",
        filtered_metadata = f"results/filtered/{output_prefix}_metadata.tsv",
        filter_log = f"logs/filter_{output_prefix}.log",
        nextclade_log = f"logs/nextclade_qc_{output_prefix}.log"
    output:
        summary_json = f"results/qc_reports/{output_prefix}_qc_summary.json",
        detailed_report = f"results/qc_reports/{output_prefix}_qc_report.html"
    params:
        segment = config["data"].get("segment", ""),
        pathogen = config["pathogen"],
        min_qc_score = config.get("qc", {}).get("nextclade", {}).get("min_qc_score", 80),
        qc_report_dir = config.get("qc", {}).get("reporting", {}).get("output_dir", "results/qc_reports")
    log:
        f"logs/generate_qc_report_{output_prefix}.log"
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p {params.qc_report_dir}
        
        # Run QC report generation script
        python scripts/generate_qc_report.py \
            --raw-sequences {input.raw_sequences} \
            --filtered-sequences {input.filtered_sequences} \
            --nextclade-json {input.nextclade_json} \
            --nextclade-passed {input.nextclade_passed} \
            --raw-metadata {input.raw_metadata} \
            --filtered-metadata {input.filtered_metadata} \
            --filter-log {input.filter_log} \
            --nextclade-log {input.nextclade_log} \
            --output-json {output.summary_json} \
            --output-html {output.detailed_report} \
            --segment {params.segment} \
            --pathogen {params.pathogen} \
            --min-qc-score {params.min_qc_score} \
            > {log} 2>&1
        """