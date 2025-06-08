"""
Rules for segment-aware sequence alignment
"""

import os

# Get filtered sequences directly (remove subsampling dependency)
def get_input_sequences(wildcards):
    return f"results/filtered/{output_prefix}_filtered.fasta"

def get_input_metadata(wildcards):
    return f"results/filtered/{output_prefix}_metadata.tsv"

def get_alignment_input(wildcards):
    """Get alignment input based on masking configuration and segment mode"""
    if segment_mode == "single":
        if config.get("mask", {}).get("sites") or config.get("mask", {}).get("from_beginning") or config.get("mask", {}).get("from_end"):
            return f"results/masked/{output_prefix}_masked.fasta"
        else:
            return f"results/aligned/{output_prefix}_aligned.fasta"
    else:
        # Multi-segment mode
        if config.get("mask", {}).get("sites") or config.get("mask", {}).get("from_beginning") or config.get("mask", {}).get("from_end"):
            return f"results/segments/{wildcards.segment}/masked/{output_prefix}_{wildcards.segment}_masked.fasta"
        else:
            return f"results/segments/{wildcards.segment}/aligned/{output_prefix}_{wildcards.segment}_aligned.fasta"

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
    threads:
        config["resources"].get("align", {}).get("threads", 4)
    log:
        "logs/align_{prefix}_{segment}.log"    
    benchmark:
        "benchmarks/align_{prefix}_{segment}.txt"
    resources:
        mem_mb = config["resources"].get("align", {}).get("mem_mb", 4000)
    shell:
        """
        python scripts/run_alignment.py --sequences {input.sequences} --output {output.alignment} --log {log} --segment {wildcards.segment} --method {params.method} --threads {threads} --mafft-options "{params.mafft_options}" --reference {params.reference}
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