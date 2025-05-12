"""
Rules for input and output validation
"""

# Validate input sequences
rule validate_fasta:
    input:
        sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta"
    output:
        validation = f"logs/validation/{output_prefix}_sequences_valid.txt"
    log:
        f"logs/validate_fasta_{output_prefix}.log"
    benchmark:
        f"benchmarks/validate_fasta_{output_prefix}.txt"
    resources:
        mem_mb = 1000
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.validation})
        
        # Validate FASTA file
        python scripts/validate_fasta.py \
            --input {input.sequences} \
            --output {output.validation} \
            > {log} 2>&1
        """

# Validate metadata
rule validate_metadata:
    input:
        metadata = f"data/metadata/raw/{output_prefix}_metadata.tsv",
        schema = "config/metadata_schema.json"
    output:
        validation = f"logs/validation/{output_prefix}_metadata_valid.txt"
    log:
        f"logs/validate_metadata_{output_prefix}.log"
    benchmark:
        f"benchmarks/validate_metadata_{output_prefix}.txt"
    resources:
        mem_mb = 1000
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.validation})
        
        # Validate metadata
        python scripts/validate_metadata.py \
            --input {input.metadata} \
            --schema {input.schema} \
            --output {output.validation} \
            > {log} 2>&1
        """

# Validate Auspice JSON output
rule validate_auspice_json:
    input:
        auspice_json = f"results/auspice/{output_prefix}.json" if segment_mode == "single" else f"results/auspice/{output_prefix}_combined.json"
    output:
        validation = f"logs/validation/{output_prefix}_auspice_valid.txt"
    log:
        f"logs/validate_auspice_json_{output_prefix}.log"
    benchmark:
        f"benchmarks/validate_auspice_json_{output_prefix}.txt"
    resources:
        mem_mb = 1000
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.validation})
        
        # Validate Auspice JSON
        python scripts/validate_auspice_json.py \
            --input {input.auspice_json} \
            --output {output.validation} \
            > {log} 2>&1
        """

# Add validation as an optional target
rule validate_all:
    input:
        fasta = f"logs/validation/{output_prefix}_sequences_valid.txt",
        metadata = f"logs/validation/{output_prefix}_metadata_valid.txt",
        auspice = f"logs/validation/{output_prefix}_auspice_valid.txt"
    output:
        touch(f"logs/validation/{output_prefix}_all_valid.txt")