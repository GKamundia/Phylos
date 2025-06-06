"""
Enhanced alignment and phylogenetic analysis rules with advanced capabilities
Integrates multiple alignment algorithms and temporal analysis
"""

# Enhanced alignment rule with multiple algorithm support
rule advanced_align:
    input:
        sequences = get_input_sequences
    output:
        alignment = f"results/aligned/{output_prefix}_aligned.fasta",
        stats = f"results/aligned/{output_prefix}_alignment_stats.json"
    params:
        reference = lambda wildcards: config["data"].get("reference_sequence", ""),
        method = config.get("align", {}).get("method", "auto"),
        threads = config["resources"].get("align", {}).get("threads", 4)
    threads: config["resources"].get("align", {}).get("threads", 4)
    log:
        f"logs/advanced_align_{output_prefix}.log"
    benchmark:
        f"benchmarks/advanced_align_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("align", {}).get("mem_mb", 4000)
    shell:
        """
        python scripts/advanced_alignment.py \
            {input.sequences} \
            {output.alignment} \
            --method {params.method} \
            --threads {params.threads} \
            --stats-output {output.stats} \
            {params.reference:+--reference {params.reference}} \
            > {log} 2>&1
        """

# Segmented alignment for multi-segment analysis
rule segment_align:
    input:
        sequences = get_input_sequences,
        metadata = get_input_metadata    output:
        alignment_dir = directory("results/segments/aligned"),
        stats = "results/segments/alignment_stats.json"
    params:
        segments = lambda wildcards: config.get("segments", ["L", "M", "S"]),
        method = config.get("align", {}).get("method", "auto"),
        threads = config["resources"].get("align", {}).get("threads", 4)
    threads: config["resources"].get("align", {}).get("threads", 4)
    log:
        "logs/segment_align.log"
    benchmark:
        f"benchmarks/segment_align_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("align", {}).get("mem_mb", 6000)
    shell:
        """
        mkdir -p {output.alignment_dir}
        
        python scripts/advanced_alignment.py \
            {input.sequences} \
            {output.alignment_dir}/combined_aligned.fasta \
            --method {params.method} \
            --threads {params.threads} \
            --metadata {input.metadata} \
            --segment-mode \
            --segments {params.segments} \
            --stats-output {output.stats} \
            > {log} 2>&1
        """

# Enhanced tree building with multiple algorithms
rule advanced_tree:
    input:
        alignment = f"results/aligned/{output_prefix}_aligned.fasta",
        metadata = get_input_metadata
    output:
        tree = f"results/tree/{output_prefix}.nwk",
        stats = f"results/tree/{output_prefix}_tree_stats.json"
    params:
        method = config.get("tree", {}).get("method", "auto"),
        model = config.get("tree", {}).get("model", "auto"),
        threads = config["resources"].get("tree", {}).get("threads", 4)
    threads: config["resources"].get("tree", {}).get("threads", 4)
    log:
        f"logs/advanced_tree_{output_prefix}.log"
    benchmark:
        f"benchmarks/advanced_tree_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("tree", {}).get("mem_mb", 8000),
        runtime = config["resources"].get("tree", {}).get("runtime", 120)
    shell:
        """
        python scripts/advanced_phylogenetics.py \
            {input.alignment} \
            {output.tree} \
            --method {params.method} \
            --model {params.model} \
            --threads {params.threads} \
            --metadata {input.metadata} \
            --stats-output {output.stats} \
            > {log} 2>&1
        """

# Temporal analysis and molecular clock estimation
rule temporal_analysis:
    input:
        tree = f"results/tree/{output_prefix}.nwk",
        metadata = get_input_metadata
    output:
        temporal_dir = directory(f"results/temporal/{output_prefix}"),
        summary = f"results/temporal/{output_prefix}/temporal_analysis_summary.txt",
        results = f"results/temporal/{output_prefix}/temporal_analysis_results.json"
    params:
        threads = config["resources"].get("temporal", {}).get("threads", 2)
    threads: config["resources"].get("temporal", {}).get("threads", 2)
    log:
        f"logs/temporal_analysis_{output_prefix}.log"
    benchmark:
        f"benchmarks/temporal_analysis_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("temporal", {}).get("mem_mb", 4000)
    shell:
        """
        mkdir -p {output.temporal_dir}
        
        python scripts/temporal_analysis.py \
            {input.metadata} \
            --tree {input.tree} \
            --output-dir {output.temporal_dir} \
            > {log} 2>&1
        """

# Molecular clock estimation with TreeTime
rule molecular_clock:
    input:
        tree = f"results/tree/{output_prefix}.nwk",
        metadata = get_input_metadata
    output:
        clock_dir = directory(f"results/molecular_clock/{output_prefix}"),
        timetree = f"results/molecular_clock/{output_prefix}/timetree.nexus"
    params:
        clock_filter = config.get("temporal", {}).get("clock_filter", 3),
        threads = config["resources"].get("temporal", {}).get("threads", 2)
    threads: config["resources"].get("temporal", {}).get("threads", 2)
    log:
        f"logs/molecular_clock_{output_prefix}.log"
    benchmark:
        f"benchmarks/molecular_clock_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("temporal", {}).get("mem_mb", 4000)
    shell:
        """
        mkdir -p {output.clock_dir}
        
        python scripts/advanced_phylogenetics.py \
            {input.tree} \
            {output.timetree} \
            --metadata {input.metadata} \
            --temporal-analysis \
            --temporal-output-dir {output.clock_dir} \
            > {log} 2>&1
        """

# Segment-specific phylogenetic analysis
rule segment_phylogenetics:
    input:
        alignment_dir = f"results/segments/aligned",
        metadata = get_input_metadata
    output:
        tree_dir = directory(f"results/segments/trees"),
        temporal_dir = directory(f"results/segments/temporal"),
        summary = f"results/segments/phylogenetic_summary.json"
    params:
        segments = lambda wildcards: config.get("segments", ["L", "M", "S"]),
        method = config.get("tree", {}).get("method", "auto"),
        threads = config["resources"].get("tree", {}).get("threads", 4)
    threads: config["resources"].get("tree", {}).get("threads", 4)
    log:
        f"logs/segment_phylogenetics_{output_prefix}.log"
    benchmark:
        f"benchmarks/segment_phylogenetics_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("tree", {}).get("mem_mb", 8000)
    run:
        import os
        import json
        
        # Create output directories
        os.makedirs(output.tree_dir, exist_ok=True)
        os.makedirs(output.temporal_dir, exist_ok=True)
        
        segment_results = {}
        
        for segment in params.segments:
            segment_alignment = os.path.join(input.alignment_dir, f"{segment}_aligned.fasta")
            
            if os.path.exists(segment_alignment):
                # Build tree for segment
                segment_tree = os.path.join(output.tree_dir, f"{segment}.nwk")
                segment_stats = os.path.join(output.tree_dir, f"{segment}_stats.json")
                
                shell(f"""
                python scripts/advanced_phylogenetics.py \
                    {segment_alignment} \
                    {segment_tree} \
                    --method {params.method} \
                    --threads {params.threads} \
                    --metadata {input.metadata} \
                    --stats-output {segment_stats} \
                    >> {log} 2>&1
                """)
                
                # Temporal analysis for segment
                segment_temporal_dir = os.path.join(output.temporal_dir, segment)
                
                shell(f"""
                mkdir -p {segment_temporal_dir}
                
                python scripts/temporal_analysis.py \
                    {input.metadata} \
                    --tree {segment_tree} \
                    --output-dir {segment_temporal_dir} \
                    >> {log} 2>&1
                """)
                
                segment_results[segment] = {
                    "tree": segment_tree,
                    "stats": segment_stats,
                    "temporal_dir": segment_temporal_dir
                }
        
        # Save summary
        with open(output.summary, 'w') as f:
            json.dump(segment_results, f, indent=2)

# Quality assessment of phylogenetic analysis
rule phylogenetic_qc:
    input:
        tree = f"results/tree/{output_prefix}.nwk",
        alignment = f"results/aligned/{output_prefix}_aligned.fasta",
        tree_stats = f"results/tree/{output_prefix}_tree_stats.json",
        temporal_results = f"results/temporal/{output_prefix}/temporal_analysis_results.json"
    output:
        qc_report = f"results/qc_reports/{output_prefix}_phylogenetic_qc.json",
        qc_summary = f"results/qc_reports/{output_prefix}_phylogenetic_qc.txt"
    log:
        f"logs/phylogenetic_qc_{output_prefix}.log"
    benchmark:
        f"benchmarks/phylogenetic_qc_{output_prefix}.txt"
    run:
        import json
        from datetime import datetime
        from Bio import Phylo, AlignIO
        
        # Initialize QC results
        qc_results = {
            "timestamp": datetime.now().isoformat(),
            "tree_qc": {},
            "alignment_qc": {},
            "temporal_qc": {},
            "overall_quality": "unknown"
        }
        
        try:
            # Tree QC
            tree = Phylo.read(input.tree, "newick")
            terminals = tree.get_terminals()
            
            qc_results["tree_qc"] = {
                "num_terminals": len(terminals),
                "tree_exists": True,
                "tree_format_valid": True
            }
            
            # Load tree statistics
            with open(input.tree_stats, 'r') as f:
                tree_stats = json.load(f)
                
            qc_results["tree_qc"].update(tree_stats.get("tree_stats", {}))
            
            # Alignment QC
            alignment = AlignIO.read(input.alignment, "fasta")
            qc_results["alignment_qc"] = {
                "num_sequences": len(alignment),
                "alignment_length": alignment.get_alignment_length(),
                "alignment_exists": True
            }
            
            # Temporal QC
            with open(input.temporal_results, 'r') as f:
                temporal_data = json.load(f)
                
            qc_results["temporal_qc"] = temporal_data.get("basic_temporal_stats", {})
            
            # Overall quality assessment
            issues = []
            
            if qc_results["tree_qc"]["num_terminals"] < 3:
                issues.append("Too few sequences for robust phylogenetic analysis")
                
            if qc_results["alignment_qc"]["num_sequences"] != qc_results["tree_qc"]["num_terminals"]:
                issues.append("Mismatch between alignment and tree sequence counts")
                
            if qc_results["temporal_qc"].get("date_range_years", 0) < 1:
                issues.append("Limited temporal range")
                
            if not issues:
                qc_results["overall_quality"] = "good"
            elif len(issues) <= 2:
                qc_results["overall_quality"] = "acceptable"
            else:
                qc_results["overall_quality"] = "poor"
                
            qc_results["issues"] = issues
            
        except Exception as e:
            qc_results["error"] = str(e)
            qc_results["overall_quality"] = "failed"
            
        # Save QC results
        with open(output.qc_report, 'w') as f:
            json.dump(qc_results, f, indent=2, default=str)
            
        # Create summary report
        with open(output.qc_summary, 'w') as f:
            f.write("PHYLOGENETIC ANALYSIS QUALITY CONTROL REPORT\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Overall Quality: {qc_results['overall_quality'].upper()}\n\n")
            
            if "tree_qc" in qc_results:
                f.write("Tree Analysis:\n")
                f.write(f"  - Number of taxa: {qc_results['tree_qc'].get('num_terminals', 'N/A')}\n")
                f.write(f"  - Tree depth: {qc_results['tree_qc'].get('tree_depth', 'N/A')}\n\n")
                
            if "alignment_qc" in qc_results:
                f.write("Alignment Analysis:\n")
                f.write(f"  - Number of sequences: {qc_results['alignment_qc'].get('num_sequences', 'N/A')}\n")
                f.write(f"  - Alignment length: {qc_results['alignment_qc'].get('alignment_length', 'N/A')}\n\n")
                
            if "temporal_qc" in qc_results:
                f.write("Temporal Analysis:\n")
                f.write(f"  - Date range: {qc_results['temporal_qc'].get('date_range_years', 'N/A'):.1f} years\n")
                f.write(f"  - Total sequences: {qc_results['temporal_qc'].get('total_sequences', 'N/A')}\n\n")
                
            if qc_results.get("issues"):
                f.write("Issues Identified:\n")
                for issue in qc_results["issues"]:
                    f.write(f"  - {issue}\n")
                f.write("\n")
                
            f.write(f"Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Rule to combine all advanced phylogenetic outputs
rule advanced_phylogenetic_summary:
    input:
        tree = f"results/tree/{output_prefix}.nwk",
        tree_stats = f"results/tree/{output_prefix}_tree_stats.json",
        temporal_summary = f"results/temporal/{output_prefix}/temporal_analysis_summary.txt",
        qc_report = f"results/qc_reports/{output_prefix}_phylogenetic_qc.json"
    output:
        summary = f"results/advanced_phylogenetic_summary_{output_prefix}.json",
        report = f"results/advanced_phylogenetic_report_{output_prefix}.txt"
    log:
        f"logs/advanced_phylogenetic_summary_{output_prefix}.log"
    run:
        import json
        from datetime import datetime
        
        # Combine all results
        summary_data = {
            "analysis_timestamp": datetime.now().isoformat(),
            "pathogen": config.get("pathogen", "unknown"),
            "output_prefix": output_prefix,
            "files_generated": {
                "tree": input.tree,
                "tree_statistics": input.tree_stats,
                "temporal_summary": input.temporal_summary,
                "qc_report": input.qc_report
            }
        }
        
        # Load and include key statistics
        try:
            with open(input.tree_stats, 'r') as f:
                tree_data = json.load(f)
                summary_data["tree_summary"] = tree_data.get("tree_stats", {})
                
            with open(input.qc_report, 'r') as f:
                qc_data = json.load(f)
                summary_data["quality_assessment"] = qc_data.get("overall_quality", "unknown")
                summary_data["issues"] = qc_data.get("issues", [])
                
        except Exception as e:
            summary_data["error"] = f"Failed to load statistics: {e}"
            
        # Save summary
        with open(output.summary, 'w') as f:
            json.dump(summary_data, f, indent=2, default=str)
            
        # Create human-readable report
        with open(output.report, 'w') as f:
            f.write("ADVANCED PHYLOGENETIC ANALYSIS REPORT\n")
            f.write("=" * 45 + "\n\n")
            f.write(f"Pathogen: {summary_data.get('pathogen', 'Unknown')}\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Quality Assessment: {summary_data.get('quality_assessment', 'Unknown').upper()}\n\n")
            
            if summary_data.get("tree_summary"):
                tree_stats = summary_data["tree_summary"]
                f.write("PHYLOGENETIC TREE STATISTICS\n")
                f.write("-" * 30 + "\n")
                f.write(f"Number of taxa: {tree_stats.get('num_terminals', 'N/A')}\n")
                f.write(f"Tree depth: {tree_stats.get('tree_depth', 'N/A')}\n")
                f.write(f"Total branch length: {tree_stats.get('total_branch_length', 'N/A')}\n")
                f.write(f"Average branch length: {tree_stats.get('avg_branch_length', 'N/A')}\n\n")
                
            if summary_data.get("issues"):
                f.write("ISSUES IDENTIFIED\n")
                f.write("-" * 18 + "\n")
                for issue in summary_data["issues"]:
                    f.write(f"• {issue}\n")
                f.write("\n")
                
            f.write("OUTPUT FILES\n")
            f.write("-" * 12 + "\n")
            for file_type, file_path in summary_data.get("files_generated", {}).items():
                f.write(f"{file_type}: {file_path}\n")
