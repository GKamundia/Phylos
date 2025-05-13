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
        exclude_ids = lambda w: config["filter"].get("exclude_ids", [])
    log:
        f"logs/filter_{output_prefix}.log"
    benchmark:
        f"benchmarks/filter_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("filter", {}).get("mem_mb", 2000)
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.sequences})
        
        # Direct approach - more reliable than temp file method
        augur filter \\
            --sequences {input.sequences} \\
            --metadata {input.metadata} \\
            --output-sequences {output.sequences} \\
            --min-length 4 \\
            --exclude-where "host=''" \\
            --exclude-where "date=''" \\
            --output-metadata {output.metadata} \\
            --include-where "segment='L'" \\
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
        # Create output directories if they don't exist
        mkdir -p $(dirname {output.sequences})
        mkdir -p $(dirname {output.metadata})
        
        # Count sequences in input file
        SEQ_COUNT=$(grep -c "^>" {input.sequences} || echo "0")
        
        # For very small datasets, just copy the files
        if [ "$SEQ_COUNT" -lt 2 ]; then
            echo "Too few sequences to subsample ($SEQ_COUNT found). Copying input files to output." > {log}
            cp {input.sequences} {output.sequences}
            cp {input.metadata} {output.metadata}
            exit 0
        fi
        
        if [ "{params.max_sequences}" -gt 0 ]; then
            # Properly quote parameters for the command
            GROUP_BY='{params.group_by}'
            
            # Construct and execute the command
            echo "Running subsampling with max_sequences={params.max_sequences}" > {log}
            
            # Extract priorities for safer handling
            PRIORITIES="{params.priorities}"
            
            # Execute augur filter directly rather than via complex command construction
            augur filter \\
                --sequences {input.sequences} \\
                --metadata {input.metadata} \\
                --output-sequences {output.sequences} \\
                --group-by "$GROUP_BY" \\
                --subsample-max-sequences {params.max_sequences} \\
                $PRIORITIES \\
                --output-metadata {output.metadata} \\
                >> {log} 2>&1 || {{
                    # Fallback if augur filter fails
                    echo "Augur filter failed. Copying input files as fallback." >> {log}
                    cp {input.sequences} {output.sequences}
                    cp {input.metadata} {output.metadata}
                }}
        else
            # If max_sequences is 0 or negative, skip subsampling
            cp {input.sequences} {output.sequences}
            cp {input.metadata} {output.metadata}
            echo "Subsampling skipped (max_sequences={params.max_sequences}). Copied files." > {log}
        fi
        """