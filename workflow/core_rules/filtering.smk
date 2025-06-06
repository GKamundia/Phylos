"""
Rules for filtering sequences and metadata
"""

# Update segment information in metadata
rule update_segment_info:
    input:
        metadata = f"data/metadata/{output_prefix}_metadata.tsv",
        sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta"
    output:
        metadata = f"data/metadata/{output_prefix}_metadata_with_segments.tsv"
    log:
        f"logs/update_segment_info_{output_prefix}.log"
    shell:
        """
        python scripts/update_segment_info.py \
            --input-metadata {input.metadata} \
            --input-sequences {input.sequences} \
            --output-metadata {output.metadata} \
            > {log} 2>&1
        """

# Filter sequences based on quality criteria and metadata
rule filter:
    input:
        sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta",
        metadata = f"data/metadata/{output_prefix}_metadata_with_segments.tsv"
    output:
        sequences = f"results/filtered/{output_prefix}_filtered.fasta",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv"
    params:
        min_length = lambda w: config["filter"].get("min_length", {}).get(
            config["data"].get("segment", "default"), 
            config["filter"].get("min_length", 0) if isinstance(config["filter"].get("min_length"), int) else 0
        ),
        max_n = config["filter"].get("max_n", 100),
        exclude_ids = lambda w: config["filter"].get("exclude_ids", [])
    log:
        f"logs/filter_{output_prefix}.log"
    benchmark:
        f"benchmarks/filter_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("filter", {}).get("mem_mb", 2000)
    shell:
        """
        python scripts/run_filter.py --sequences {input.sequences} --metadata {input.metadata} --output-sequences {output.sequences} --output-metadata {output.metadata} --log {log}
        """

# Subsample sequences (optional, based on configuration)
rule subsample:
    input:
        sequences = f"results/filtered/{output_prefix}_filtered.fasta",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv"
    output:
        sequences = f"results/subsampled/{output_prefix}_subsampled.fasta",
        metadata = f"results/subsampled/{output_prefix}_metadata.tsv"
    params:
        group_by = config["subsample"].get("group_by", "country year"),
        max_sequences = config["subsample"].get("max_sequences", 1000),
        priorities = lambda w: " ".join([
            f"--priority {priority_type} {field}"
            for priority_type, field in {
                "recency": "num_date",
                **{k: v.get("field") for k, v in config["subsample"].get("priorities", {}).items() if v.get("type") == "numeric"}
            }.items()
        ])
    log:
        f"logs/subsample_{output_prefix}.log"
    benchmark:
        f"benchmarks/subsample_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("subsample", {}).get("mem_mb", 2000)
    shell:
        """
        python scripts/run_subsample.py --input-sequences {input.sequences} --input-metadata {input.metadata} --output-sequences {output.sequences} --output-metadata {output.metadata} --log {log}
        """