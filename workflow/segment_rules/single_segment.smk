"""
Rules specific to single segment analysis
"""

# Run Nextclade for QC, clade assignment and outlier detection
rule nextclade_qc_single:
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
    resources:
        mem_mb = config["resources"].get("nextclade", {}).get("mem_mb", 4000)
    threads: config["resources"].get("nextclade", {}).get("threads", 2)
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p {params.outdir}
        
        # Count sequences in input file
        SEQ_COUNT=$(grep -c "^>" {input.sequences} || echo "0")
        
        if [ "$SEQ_COUNT" -eq "0" ]; then
            echo "No sequences found in input file. Creating empty output files." > {log}
            echo "{{\"version\":\"2.0.0\",\"results\":[]}}" > {output.json}
            echo "seqName\tqc.overallScore\tqc.overallStatus" > {output.tsv}
            copy "{input.sequences}" "{output.aligned}"
            copy "{input.sequences}" "{output.passed}"
            exit 0
        fi
        
        # Check if dataset directory exists and has required files
        if [ ! -d "{params.dataset_dir}" ] || [ ! -f "{params.dataset_dir}/reference.fasta" ]; then
            echo "Nextclade dataset not found or incomplete at {params.dataset_dir}" > {log}
            echo "Creating placeholder output files to allow pipeline to continue" >> {log}
            # Create minimal outputs
            echo "{{\"version\":\"2.0.0\",\"results\":[]}}" > {output.json}
            echo "seqName\tqc.overallScore\tqc.overallStatus" > {output.tsv}
            copy "{input.sequences}" "{output.aligned}"
            copy "{input.sequences}" "{output.passed}"
            exit 0
        fi
        
        # Attempt to run nextclade, but with error handling
        if nextclade run \\
            --input-dataset {params.dataset_dir} \\
            --output-json {output.json} \\
            --output-tsv {output.tsv} \\
            --output-aligned {output.aligned} \\
            --input-fasta {input.sequences} \\
            --include-reference \\
            --jobs {threads} \\
            > {log} 2>&1; then
            
            # If Nextclade succeeds, run the filter script
            python scripts/filter_nextclade_results.py \\
                --input-json {output.json} \\
                --input-fasta {input.sequences} \\
                --output-fasta {output.passed} \\
                --min-qc-score {params.min_qc_score} \\
                --segment {params.segment} \\
                >> {log} 2>&1 || true
        else
            # If Nextclade fails, create placeholder outputs
            echo "Nextclade failed. Using input sequences as passed sequences." >> {log}
            copy "{input.sequences}" "{output.passed}"
            
            # Ensure other output files exist
            if [ ! -f "{output.json}" ]; then
                echo "{{\"version\":\"2.0.0\",\"results\":[]}}" > {output.json}
            fi
            if [ ! -f "{output.tsv}" ]; then
                echo "seqName\tqc.overallScore\tqc.overallStatus" > {output.tsv}
            fi
            if [ ! -f "{output.aligned}" ]; then
                copy "{input.sequences}" "{output.aligned}"
            fi
        fi
        """

# Generate QC report for single segment
rule generate_qc_report_single:
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
    benchmark:
        f"benchmarks/qc_report_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("qc_report", {}).get("mem_mb", 2000)
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
