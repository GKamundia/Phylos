"""
Advanced Phylogenetics Rules for Enhanced Analysis
Provides sophisticated phylogenetic analysis capabilities including
multiple alignment methods, advanced tree building, and temporal analysis.
"""

rule advanced_align:
    """
    Enhanced sequence alignment with multiple algorithm support
    """
    input:
        sequences="results/filtered/{pathogen}_filtered.fasta"
    output:
        alignment="results/advanced_phylogenetics/{pathogen}_advanced_alignment.fasta",
        report="results/advanced_phylogenetics/{pathogen}_alignment_report.json"
    params:
        method=lambda w: config.get("advanced_phylogenetics", {}).get("alignment", {}).get("method", "mafft"),
        threads=lambda w: config.get("advanced_phylogenetics", {}).get("alignment", {}).get("threads", 4),
        additional_params=lambda w: config.get("advanced_phylogenetics", {}).get("alignment", {}).get("additional_params", "")
    resources:
        mem_mb=lambda w: config.get("resources", {}).get("alignment_memory", 8000),
        runtime=lambda w: config.get("resources", {}).get("alignment_time", 120)
    benchmark:
        "benchmarks/advanced_align_{pathogen}.txt"
    log:
        "logs/advanced_align_{pathogen}.log"
    run:
        try:
            import sys
            sys.path.append("scripts")
            from advanced_alignment import AdvancedAligner
            
            aligner = AdvancedAligner(
                method=params.method,
                threads=params.threads,
                additional_params=params.additional_params
            )
            
            aligner.align_sequences(
                input_fasta=input.sequences,
                output_fasta=output.alignment,
                report_file=output.report
            )
            
        except Exception as e:
            with open(log[0], "w") as log_file:
                log_file.write(f"Advanced alignment failed: {str(e)}\n")
            raise

rule advanced_tree:
    """
    Multi-method tree building with bootstrap support and model selection
    """
    input:
        alignment="results/advanced_phylogenetics/{pathogen}_advanced_alignment.fasta"
    output:
        tree="results/advanced_phylogenetics/{pathogen}_advanced_tree.nwk",
        bootstrap="results/advanced_phylogenetics/{pathogen}_bootstrap_trees.nwk",
        report="results/advanced_phylogenetics/{pathogen}_tree_report.json"
    params:
        method=lambda w: config.get("advanced_phylogenetics", {}).get("tree_building", {}).get("method", "iqtree"),
        model=lambda w: config.get("advanced_phylogenetics", {}).get("tree_building", {}).get("model", "GTR+G"),
        bootstrap_replicates=lambda w: config.get("advanced_phylogenetics", {}).get("tree_building", {}).get("bootstrap_replicates", 1000),
        threads=lambda w: config.get("advanced_phylogenetics", {}).get("tree_building", {}).get("threads", 4)
    resources:
        mem_mb=lambda w: config.get("resources", {}).get("tree_memory", 16000),
        runtime=lambda w: config.get("resources", {}).get("tree_time", 240)
    benchmark:
        "benchmarks/advanced_tree_{pathogen}.txt"
    log:
        "logs/advanced_tree_{pathogen}.log"
    run:
        try:
            import sys
            sys.path.append("scripts")
            from advanced_phylogenetics import AdvancedTreeBuilder
            
            tree_builder = AdvancedTreeBuilder(
                method=params.method,
                model=params.model,
                bootstrap_replicates=params.bootstrap_replicates,
                threads=params.threads
            )
            
            tree_builder.build_tree(
                alignment_file=input.alignment,
                output_tree=output.tree,
                bootstrap_trees=output.bootstrap,
                report_file=output.report
            )
            
        except Exception as e:
            with open(log[0], "w") as log_file:
                log_file.write(f"Advanced tree building failed: {str(e)}\n")
            raise

rule temporal_analysis:
    """
    Comprehensive temporal pattern analysis and molecular clock estimation
    """    input:
        tree="results/advanced_phylogenetics/{pathogen}_advanced_tree.nwk",
        metadata="results/filtered/{pathogen}_metadata.tsv"
    output:
        temporal_tree="results/advanced_phylogenetics/{pathogen}_temporal_tree.nwk",
        clock_report="results/advanced_phylogenetics/{pathogen}_clock_analysis.json",
        temporal_plots="results/advanced_phylogenetics/{pathogen}_temporal_plots.pdf"
    params:
        clock_model=lambda w: config.get("advanced_phylogenetics", {}).get("temporal_analysis", {}).get("clock_model", "strict"),
        date_inference=lambda w: config.get("advanced_phylogenetics", {}).get("temporal_analysis", {}).get("date_inference", True),
        outlier_detection=lambda w: config.get("advanced_phylogenetics", {}).get("temporal_analysis", {}).get("outlier_detection", True)
    resources:
        mem_mb=lambda w: config.get("resources", {}).get("temporal_memory", 8000),
        runtime=lambda w: config.get("resources", {}).get("temporal_time", 180)
    benchmark:
        "benchmarks/temporal_analysis_{pathogen}.txt"
    log:
        "logs/temporal_analysis_{pathogen}.log"
    run:
        try:
            import sys
            sys.path.append("scripts")
            from temporal_analysis import TemporalAnalyzer
            
            analyzer = TemporalAnalyzer(
                clock_model=params.clock_model,
                date_inference=params.date_inference,
                outlier_detection=params.outlier_detection
            )
            
            analyzer.analyze_temporal_patterns(
                tree_file=input.tree,
                metadata_file=input.metadata,
                output_tree=output.temporal_tree,
                clock_report=output.clock_report,
                plots_file=output.temporal_plots
            )
            
        except Exception as e:
            with open(log[0], "w") as log_file:
                log_file.write(f"Temporal analysis failed: {str(e)}\n")
            raise

rule phylogenetic_qc:
    """
    Quality control assessment for phylogenetic analysis
    """
    input:
        alignment="results/advanced_phylogenetics/{pathogen}_advanced_alignment.fasta",
        tree="results/advanced_phylogenetics/{pathogen}_advanced_tree.nwk",
        bootstrap="results/advanced_phylogenetics/{pathogen}_bootstrap_trees.nwk"
    output:
        qc_report="results/advanced_phylogenetics/{pathogen}_phylo_qc.json",
        qc_plots="results/advanced_phylogenetics/{pathogen}_phylo_qc_plots.pdf"
    params:
        quality_threshold=lambda w: config.get("advanced_phylogenetics", {}).get("quality_control", {}).get("alignment_quality_threshold", 0.7),
        bootstrap_threshold=lambda w: config.get("advanced_phylogenetics", {}).get("quality_control", {}).get("bootstrap_threshold", 70),
        branch_length_threshold=lambda w: config.get("advanced_phylogenetics", {}).get("quality_control", {}).get("branch_length_threshold", 0.01)
    resources:
        mem_mb=lambda w: config.get("resources", {}).get("qc_memory", 4000),
        runtime=lambda w: config.get("resources", {}).get("qc_time", 60)
    benchmark:
        "benchmarks/phylogenetic_qc_{pathogen}.txt"
    log:
        "logs/phylogenetic_qc_{pathogen}.log"
    run:
        try:
            import sys
            sys.path.append("scripts")
            from advanced_phylogenetics import PhylogeneticQC
            
            qc_analyzer = PhylogeneticQC(
                quality_threshold=params.quality_threshold,
                bootstrap_threshold=params.bootstrap_threshold,
                branch_length_threshold=params.branch_length_threshold
            )
            
            qc_analyzer.assess_quality(
                alignment_file=input.alignment,
                tree_file=input.tree,
                bootstrap_file=input.bootstrap,
                report_file=output.qc_report,
                plots_file=output.qc_plots
            )
            
        except Exception as e:
            with open(log[0], "w") as log_file:
                log_file.write(f"Phylogenetic QC failed: {str(e)}\n")
            raise

rule advanced_phylogenetic_report:
    """
    Generate comprehensive HTML and JSON reports for advanced phylogenetic analysis
    """
    input:
        alignment_report="results/advanced_phylogenetics/{pathogen}_alignment_report.json",
        tree_report="results/advanced_phylogenetics/{pathogen}_tree_report.json",
        clock_report="results/advanced_phylogenetics/{pathogen}_clock_analysis.json",
        qc_report="results/advanced_phylogenetics/{pathogen}_phylo_qc.json"
    output:
        html_report="results/advanced_phylogenetics/{pathogen}_advanced_phylo_report.html",
        json_summary="results/advanced_phylogenetics/{pathogen}_advanced_phylo_summary.json"
    params:
        report_title=lambda w: f"Advanced Phylogenetic Analysis - {config.get('pathogen_name', w.pathogen)}",
        include_plots=lambda w: config.get("advanced_phylogenetics", {}).get("reporting", {}).get("include_plots", True)
    resources:
        mem_mb=lambda w: config.get("resources", {}).get("report_memory", 2000),
        runtime=lambda w: config.get("resources", {}).get("report_time", 30)
    benchmark:
        "benchmarks/advanced_phylogenetic_report_{pathogen}.txt"
    log:
        "logs/advanced_phylogenetic_report_{pathogen}.log"
    run:
        try:
            import sys
            sys.path.append("scripts")
            from run_advanced_phylogenetics import generate_comprehensive_report
            
            generate_comprehensive_report(
                alignment_report=input.alignment_report,
                tree_report=input.tree_report,
                clock_report=input.clock_report,
                qc_report=input.qc_report,
                html_output=output.html_report,
                json_output=output.json_summary,
                title=params.report_title,
                include_plots=params.include_plots
            )
            
        except Exception as e:
            with open(log[0], "w") as log_file:
                log_file.write(f"Report generation failed: {str(e)}\n")
            raise

# Target rule for complete advanced phylogenetic analysis
rule run_advanced_phylogenetics:
    """
    Complete advanced phylogenetic analysis pipeline
    """
    input:
        "results/advanced_phylogenetics/{pathogen}_advanced_phylo_report.html",
        "results/advanced_phylogenetics/{pathogen}_advanced_phylo_summary.json"
    output:
        touch("results/advanced_phylogenetics/.{pathogen}_analysis_complete")
    message:
        "Advanced phylogenetic analysis completed for {wildcards.pathogen}"