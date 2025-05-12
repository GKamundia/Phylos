"""
Common rules shared between single and multi-segment workflows
"""

# Reconstruct ancestral traits (e.g., country)
rule traits:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/tree/{output_prefix}_{wildcards.segment}_refined.nwk",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/filtered/{output_prefix}_{wildcards.segment}_metadata.tsv"
    output:
        node_data = f"results/node_data/{output_prefix}_traits.json" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/node_data/{output_prefix}_{wildcards.segment}_traits.json"
    params:
        columns = lambda w: ",".join(config.get("traits", {}).get("columns", ["country"]))
    log:
        f"logs/traits_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/traits_{output_prefix}_{wildcards.segment}.log"
    benchmark:
        f"benchmarks/traits_{output_prefix}.txt" if segment_mode == "single" else lambda wildcards: f"benchmarks/traits_{output_prefix}_{wildcards.segment}.txt"
    resources:
        mem_mb = config["resources"].get("traits", {}).get("mem_mb", 2000)
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            > {log} 2>&1
        """

# Export to auspice
rule export:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/tree/{output_prefix}_{wildcards.segment}_refined.nwk",
        metadata = f"results/filtered/{output_prefix}_metadata.tsv" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/filtered/{output_prefix}_{wildcards.segment}_metadata.tsv",
        branch_lengths = f"results/node_data/{output_prefix}_branch_lengths.json" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/node_data/{output_prefix}_{wildcards.segment}_branch_lengths.json",
        traits = f"results/node_data/{output_prefix}_traits.json" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/node_data/{output_prefix}_{wildcards.segment}_traits.json"
    output:
        auspice_json = f"results/auspice/{output_prefix}.json" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/auspice/{output_prefix}_{wildcards.segment}.json"
    params:
        auspice_config = pathogen_config.get("auspice_config_path", "config/auspice_config.json"),
        lat_longs = pathogen_config.get("lat_longs_path", "config/lat_longs.tsv"),
        title = lambda wildcards: f"{config['pathogen_name']} Genomic Surveillance" if segment_mode == "single" else f"{config['pathogen_name']} {wildcards.segment} Segment Surveillance"
    log:
        f"logs/export_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/export_{output_prefix}_{wildcards.segment}.log"
    benchmark:
        f"benchmarks/export_{output_prefix}.txt" if segment_mode == "single" else lambda wildcards: f"benchmarks/export_{output_prefix}_{wildcards.segment}.txt"
    resources:
        mem_mb = config["resources"].get("export", {}).get("mem_mb", 2000)
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} \
            --auspice-config {params.auspice_config} \
            --lat-longs {params.lat_longs} \
            --title "{params.title}" \
            --output {output.auspice_json} \
            > {log} 2>&1
        """

# Run Nextclade for sequence QC, clade assignment, and outlier detection
rule nextclade_qc:
    input:
        sequences = f"results/filtered/{output_prefix}_filtered.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/filtered/{output_prefix}_{wildcards.segment}_filtered.fasta"
    output:
        json = f"results/nextclade/{output_prefix}_nextclade.json" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/nextclade/{output_prefix}_{wildcards.segment}_nextclade.json",
        tsv = f"results/nextclade/{output_prefix}_nextclade.tsv" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/nextclade/{output_prefix}_{wildcards.segment}_nextclade.tsv",
        aligned = f"results/nextclade/{output_prefix}_aligned.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/nextclade/{output_prefix}_{wildcards.segment}_aligned.fasta",
        passed = f"results/nextclade/{output_prefix}_passed.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/nextclade/{output_prefix}_{wildcards.segment}_passed.fasta"
    params:
        segment = config["data"].get("segment", "") if segment_mode == "single" else lambda wildcards: wildcards.segment,
        dataset_dir = lambda w: f"nextclade/datasets/{config['pathogen']}/",
        min_qc_score = config.get("qc", {}).get("nextclade", {}).get("min_qc_score", 80),
        outdir = "results/nextclade" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/nextclade"
    log:
        f"logs/nextclade_qc_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/nextclade_qc_{output_prefix}_{wildcards.segment}.log"
    benchmark:
        f"benchmarks/nextclade_qc_{output_prefix}.txt" if segment_mode == "single" else lambda wildcards: f"benchmarks/nextclade_qc_{output_prefix}_{wildcards.segment}.txt"
    resources:
        mem_mb = config["resources"].get("nextclade", {}).get("mem_mb", 4000)
    threads: config["resources"].get("nextclade", {}).get("threads", 2)
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
            --jobs {threads} \
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

# Generate QC report
rule generate_qc_report:
    input:
        raw_sequences = f"data/sequences/raw/{output_prefix}_sequences.fasta" if segment_mode == "single" else lambda wildcards: f"data/sequences/raw/{output_prefix}_sequences.fasta",
        filtered_sequences = f"results/filtered/{output_prefix}_filtered.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/filtered/{output_prefix}_{wildcards.segment}_filtered.fasta",
        nextclade_json = f"results/nextclade/{output_prefix}_nextclade.json" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/nextclade/{output_prefix}_{wildcards.segment}_nextclade.json",
        nextclade_passed = f"results/nextclade/{output_prefix}_passed.fasta" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/nextclade/{output_prefix}_{wildcards.segment}_passed.fasta",
        raw_metadata = f"data/metadata/raw/{output_prefix}_metadata.tsv" if segment_mode == "single" else lambda wildcards: f"data/metadata/raw/{output_prefix}_metadata.tsv",
        filtered_metadata = f"results/filtered/{output_prefix}_metadata.tsv" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/filtered/{output_prefix}_{wildcards.segment}_metadata.tsv",
        filter_log = f"logs/filter_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/filter_{output_prefix}_{wildcards.segment}.log",
        nextclade_log = f"logs/nextclade_qc_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/nextclade_qc_{output_prefix}_{wildcards.segment}.log"
    output:
        summary_json = f"results/qc_reports/{output_prefix}_qc_summary.json" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/qc_reports/{output_prefix}_{wildcards.segment}_qc_summary.json",
        detailed_report = f"results/qc_reports/{output_prefix}_qc_report.html" if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/qc_reports/{output_prefix}_{wildcards.segment}_qc_report.html"
    params:
        segment = config["data"].get("segment", "") if segment_mode == "single" else lambda wildcards: wildcards.segment,
        pathogen = config["pathogen"],
        min_qc_score = config.get("qc", {}).get("nextclade", {}).get("min_qc_score", 80),
        qc_report_dir = config.get("qc", {}).get("reporting", {}).get("output_dir", "results/qc_reports") if segment_mode == "single" else lambda wildcards: f"results/segments/{wildcards.segment}/qc_reports"
    log:
        f"logs/generate_qc_report_{output_prefix}.log" if segment_mode == "single" else lambda wildcards: f"logs/generate_qc_report_{output_prefix}_{wildcards.segment}.log"
    benchmark:
        f"benchmarks/qc_report_{output_prefix}.txt" if segment_mode == "single" else lambda wildcards: f"benchmarks/qc_report_{output_prefix}_{wildcards.segment}.txt"
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