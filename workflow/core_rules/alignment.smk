"""
Rules for segment-aware sequence alignment
"""

# Get filtered sequences directly (remove subsampling dependency)
def get_input_sequences(wildcards):
    return f"results/filtered/{output_prefix}_filtered.fasta"

def get_input_metadata(wildcards):
    return f"results/filtered/{output_prefix}_metadata.tsv"

# Split sequences by segment first
checkpoint split_by_segment:
    input:
        sequences = get_input_sequences,
        metadata = get_input_metadata
    output:
        sequences = expand("results/segments/{segment}/raw/{output_prefix}_{segment}_sequences.fasta", segment=["L", "M", "S"]),
        metadata = expand("results/segments/{segment}/filtered/{output_prefix}_{segment}_metadata.tsv", segment=["L", "M", "S"])
    params:
        segments = ["L", "M", "S"]
    log:
        f"logs/split_segments_{output_prefix}.log"
    benchmark:
        f"benchmarks/split_segments_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("split_segments", {}).get("mem_mb", 2000)
    script:
        "../../scripts/split_by_segment.py"

# Align each segment separately
rule align_segment:
    input:
        sequences = "results/segments/{segment}/raw/{prefix}_{segment}_sequences.fasta"
    output:
        alignment = "results/segments/{segment}/aligned/{prefix}_{segment}_aligned.fasta"
    params:
        reference = lambda wildcards: get_segment_reference(wildcards.segment),
        method = config.get("align", {}).get("method", "mafft"),
        mafft_options = config.get("align", {}).get("mafft_options", "--nomemsave --auto")
    threads: config["resources"].get("align", {}).get("threads", 4)
    log:
        "logs/align_{prefix}_{segment}.log"
    benchmark:
        "benchmarks/align_{prefix}_{segment}.txt"
    resources:
        mem_mb = config["resources"].get("align", {}).get("mem_mb", 4000)
    shell:
        """
        mkdir -p $(dirname {output.alignment})
        
        # Check if we have sequences to align
        if [ ! -s {input.sequences} ]; then
            echo "No sequences found for segment {wildcards.segment}" > {log}
            touch {output.alignment}
            exit 0
        fi
        
        # Count sequences
        seq_count=$(grep -c "^>" {input.sequences} || echo "0")
        echo "Aligning $seq_count sequences for segment {wildcards.segment}" > {log}
        
        if [ "$seq_count" -eq 0 ]; then
            echo "No sequences to align for segment {wildcards.segment}" >> {log}
            touch {output.alignment}
            exit 0
        elif [ "$seq_count" -eq 1 ]; then
            echo "Only one sequence found, copying as alignment" >> {log}
            cp {input.sequences} {output.alignment}
        else
            if [ "{params.method}" = "mafft" ]; then
                if [ -n "{params.reference}" ] && [ -f "{params.reference}" ]; then
                    echo "Using reference sequence: {params.reference}" >> {log}
                    augur align --sequences {input.sequences} --output {output.alignment} --reference-sequence {params.reference} --nthreads {threads} --method mafft --mafft-options "{params.mafft_options}" >> {log} 2>&1
                else
                    echo "No reference sequence, aligning without reference" >> {log}
                    augur align --sequences {input.sequences} --output {output.alignment} --nthreads {threads} --method mafft --mafft-options "{params.mafft_options}" >> {log} 2>&1
                fi
            else
                if [ -n "{params.reference}" ] && [ -f "{params.reference}" ]; then
                    augur align --sequences {input.sequences} --output {output.alignment} --reference-sequence {params.reference} --nthreads {threads} >> {log} 2>&1
                else
                    augur align --sequences {input.sequences} --output {output.alignment} --nthreads {threads} >> {log} 2>&1
                fi
            fi
        fi
        """

# Helper function to get segment-specific reference sequences
def get_segment_reference(segment):
    """Get NCBI RefSeq reference sequence for specific segment"""
    # Use NCBI RefSeq references extracted from raw data
    ref_path = f"data/references/reference_{segment.lower()}.fasta"
    
    # Check if reference file exists
    if os.path.exists(ref_path):
        return ref_path
    else:
        print(f"Warning: Reference sequence not found for segment {segment} at {ref_path}")
        return ""

# Mask sites in segment alignment (optional)
rule mask_segment:
    input:
        alignment = "results/segments/{segment}/aligned/{prefix}_{segment}_aligned.fasta"
    output:
        alignment = "results/segments/{segment}/masked/{prefix}_{segment}_masked.fasta"
    params:
        mask_sites = lambda w: " ".join(
            f"--mask-sites {site}" for site in config.get("mask", {}).get("sites", {}).get(w.segment, [])
        ),
        mask_from = lambda w: f"--mask-from-beginning {config['mask']['from_beginning'][w.segment]}" 
                            if config.get("mask", {}).get("from_beginning", {}).get(w.segment) else "",
        mask_to = lambda w: f"--mask-from-end {config['mask']['from_end'][w.segment]}" 
                          if config.get("mask", {}).get("from_end", {}).get(w.segment) else ""
    log:
        "logs/mask_{prefix}_{segment}.log"
    benchmark:
        "benchmarks/mask_{prefix}_{segment}.txt"
    resources:
        mem_mb = config["resources"].get("mask", {}).get("mem_mb", 2000)
    shell:
        """
        if [ ! -s {input.alignment} ]; then
            echo "No alignment to mask for segment {wildcards.segment}" > {log}
            touch {output.alignment}
        elif [ -n "{params.mask_sites}" ] || [ -n "{params.mask_from}" ] || [ -n "{params.mask_to}" ]; then
            augur mask --sequences {input.alignment} --output {output.alignment} {params.mask_sites} {params.mask_from} {params.mask_to} > {log} 2>&1
        else
            cp {input.alignment} {output.alignment}
        fi
        """