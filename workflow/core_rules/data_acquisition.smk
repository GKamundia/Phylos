"""
Rules for data acquisition and initial metadata preparation
"""

# Download sequence data from source
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
    shell:
        """
        python scripts/download_ncbi_virus_exact.py \
            --email "{params.email}" \
            --output-fasta {output.sequences} \
            --output-metadata {output.metadata} \
            > {log} 2>&1
        """

# Prepare and validate metadata
checkpoint prepare_metadata:    
    input:
        metadata = f"data/metadata/raw/{output_prefix}_metadata.tsv",
        schema = "config/metadata_schema.json",
        lat_longs = config.get("lat_longs_path", "config/lat_longs.tsv")
    output:
        metadata = f"data/metadata/{output_prefix}_metadata.tsv",
        report = f"data/metadata/{output_prefix}_validation_report.json"
    params:
        strict_flag = "--strict" if config.get("workflow", {}).get("strict_metadata", False) else ""
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
            {params.strict_flag} \
            > {log} 2>&1
        """