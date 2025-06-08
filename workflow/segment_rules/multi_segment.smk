"""
Rules specific to handling multiple genome segments
"""

# Use checkpoint to determine how to split segments
checkpoint split_by_segment:
    input:
        sequences = "results/filtered/{prefix}_filtered.fasta",
        metadata = "results/filtered/{prefix}_metadata.tsv"
    output:
        sequences = expand("results/segments/{segment}/raw/{{prefix}}_{segment}_sequences.fasta", segment=segments),
        metadata = expand("results/segments/{segment}/filtered/{{prefix}}_{segment}_metadata.tsv", segment=segments),
        done = "results/segments/split_by_segment_{prefix}.done"
    params:
        segments = segments
    log:
        "logs/split_segments_{prefix}.log"
    benchmark:
        "benchmarks/split_segments_{prefix}.txt"
    resources:
        mem_mb = config["resources"].get("split_segments", {}).get("mem_mb", 2000)
    script:
        "../../scripts/split_by_segment.py"

# Generate combined visualization for all segments
rule combine_segments:
    input:
        segment_jsons = lambda wildcards: expand("results/segments/{segment}/auspice/{{prefix}}_{segment}.json", segment=get_available_segments(wildcards))
    output:
        combined_json = "results/auspice/{prefix}_combined.json"
    params:
        title = f"{config['pathogen_name']} Multi-Segment Analysis",
        segments = segments
    log:
        "logs/combine_segments_{prefix}.log"
    benchmark:
        "benchmarks/combine_segments_{prefix}.txt"
    resources:
        mem_mb = config["resources"].get("combine", {}).get("mem_mb", 2000)
    script:
        "../../scripts/combine_segment_jsons.py"

# Generate combined QC report for all segments
rule combined_qc_report:
    input:
        segment_qc = lambda wildcards: expand("results/segments/{segment}/qc_reports/{{prefix}}_{segment}_qc_summary.json", segment=get_available_segments(wildcards))
    output:
        combined_report = "results/qc_reports/{prefix}_qc_summary.json"
    params:
        pathogen = config["pathogen"],
        segments = segments
    log:
        "logs/combined_qc_report_{prefix}.log"
    benchmark:
        "benchmarks/combined_qc_report_{prefix}.txt"
    resources:
        mem_mb = config["resources"].get("combined_qc", {}).get("mem_mb", 2000)
    script:
        "../../scripts/generate_combined_qc_report.py"