"""
Rules for filtering sequences and metadata
"""

# Filter sequences based on quality criteria and metadata
rule filter:
    input:
        sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta",
        metadata = f"data/metadata/{output_prefix}_metadata.tsv"
    output:
        sequences = f"results/filtered/{output_prefix}_filtered.fasta",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv"
    params:
        min_length = lambda w: config["filter"].get("min_length", {}).get(
            config["data"].get("segment", "default"), 
            config["filter"].get("min_length", 0) if isinstance(config["filter"].get("min_length"), int) else 0
        ),
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
        if [ "{params.max_sequences}" -gt 0 ]; then
            augur filter \
                --sequences {input.sequences} \
                --metadata {input.metadata} \
                --output {output.sequences} \
                --group-by {params.group_by} \
                --subsample-max-sequences {params.max_sequences} \
                {params.priorities} \
                --output-metadata {output.metadata} \
                > {log} 2>&1
        else
            # If max_sequences is 0 or negative, skip subsampling
            cp {input.sequences} {output.sequences}
            cp {input.metadata} {output.metadata}
        fi
        """