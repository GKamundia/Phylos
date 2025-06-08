"""
Rules for sequence alignment
"""

# Get filtered or subsampled sequences based on configuration
def get_input_sequences(wildcards):
    if config["subsample"].get("max_sequences", 0) > 0:
        return f"results/subsampled/{output_prefix}_subsampled.fasta"
    else:
        return f"results/filtered/{output_prefix}_filtered.fasta"

def get_input_metadata(wildcards):
    if config["subsample"].get("max_sequences", 0) > 0:
        return f"results/subsampled/{output_prefix}_metadata.tsv"
    else:
        # Always use filtered metadata (segment fixes are now in prepare_metadata)
        return f"results/filtered/{output_prefix}_metadata.tsv"

def get_alignment_input(wildcards):
    """Get alignment input based on masking configuration"""
    if config.get("mask", {}).get("sites") or config.get("mask", {}).get("from_beginning") or config.get("mask", {}).get("from_end"):
        return f"results/masked/{output_prefix}_masked.fasta"
    else:
        return f"results/aligned/{output_prefix}_aligned.fasta"

# Align sequences
rule align:
    input:
        sequences = get_input_sequences
    output:
        alignment = f"results/aligned/{output_prefix}_aligned.fasta"
    params:
        reference = lambda wildcards: config["data"].get("reference_sequence", ""),
        method = config.get("align", {}).get("method", "mafft"),
        mafft_options = config.get("align", {}).get("mafft_options", "--nomemsave --auto")
    threads: config["resources"].get("align", {}).get("threads", 4)
    log:
        f"logs/align_{output_prefix}.log"
    benchmark:
        f"benchmarks/align_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("align", {}).get("mem_mb", 4000)
    shell:
        """
        mkdir -p $(dirname {output.alignment})
        
        if [ "{params.method}" = "mafft" ]; then
            if [ -n "{params.reference}" ]; then
                augur align --sequences {input.sequences} --output {output.alignment} --reference-sequence {params.reference} --nthreads {threads} --method mafft --mafft-options "{params.mafft_options}" > {log} 2>&1
            else
                augur align --sequences {input.sequences} --output {output.alignment} --nthreads {threads} --method mafft --mafft-options "{params.mafft_options}" > {log} 2>&1
            fi
        else
            if [ -n "{params.reference}" ]; then
                augur align --sequences {input.sequences} --output {output.alignment} --reference-sequence {params.reference} --nthreads {threads} > {log} 2>&1
            else
                augur align --sequences {input.sequences} --output {output.alignment} --nthreads {threads} > {log} 2>&1
            fi
        fi
        """

# Mask sites in alignment (optional)
rule mask:
    input:
        alignment = f"results/aligned/{output_prefix}_aligned.fasta"
    output:
        alignment = f"results/masked/{output_prefix}_masked.fasta"
    params:
        mask_sites = lambda w: " ".join(
            f"--mask-sites {site}" for site in config.get("mask", {}).get("sites", [])
        ),
        mask_from = lambda w: f"--mask-from-beginning {config['mask']['from_beginning']}" 
                            if config.get("mask", {}).get("from_beginning") else "",
        mask_to = lambda w: f"--mask-from-end {config['mask']['from_end']}" 
                          if config.get("mask", {}).get("from_end") else ""
    log:
        f"logs/mask_{output_prefix}.log"
    benchmark:
        f"benchmarks/mask_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("mask", {}).get("mem_mb", 2000)
    shell:
        """        if [ -n "{params.mask_sites}" ] || [ -n "{params.mask_from}" ] || [ -n "{params.mask_to}" ]; then
            augur mask --sequences {input.alignment} --output {output.alignment} {params.mask_sites} {params.mask_from} {params.mask_to} > {log} 2>&1
        else
            copy "{input.alignment}" "{output.alignment}"
        fi
        """