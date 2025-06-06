"""
Advanced Phylogenetics Rules for Nextstrain Pipeline

This module provides enhanced phylogenetic analysis capabilities including:
- Advanced alignment algorithms
- Multiple tree building methods  
- Temporal analysis
- Molecular clock estimation
"""

# Enhanced alignment rule with multiple algorithm support
rule advanced_align:
    input:
        sequences = get_input_sequences
    output:
        alignment = "results/aligned/{pathogen}_aligned.fasta",
        stats = "results/aligned/{pathogen}_alignment_stats.json"
    params:
        method = lambda wildcards: config.get("advanced_phylogenetics", {}).get("alignment", {}).get("methods", ["auto"])[0],
        threads = lambda wildcards: config["resources"].get("align", {}).get("threads", 4)
    threads: lambda wildcards: config["resources"].get("align", {}).get("threads", 4)
    log:
        "logs/advanced_align_{pathogen}.log"
    benchmark:
        "benchmarks/advanced_align_{pathogen}.txt"
    resources:
        mem_mb = lambda wildcards: config["resources"].get("align", {}).get("mem_mb", 4000)
    shell:
        """
        python scripts/advanced_alignment.py \
            {input.sequences} \
            {output.alignment} \
            --method {params.method} \
            --threads {params.threads} \
            --stats-output {output.stats} \
            > {log} 2>&1
        """

# Enhanced tree building with multiple algorithms
rule advanced_tree:
    input:
        alignment = "results/aligned/{pathogen}_aligned.fasta",
        metadata = get_input_metadata
    output:
        tree = "results/tree/{pathogen}.nwk",
        stats = "results/tree/{pathogen}_tree_stats.json"
    params:
        method = lambda wildcards: config.get("advanced_phylogenetics", {}).get("tree_building", {}).get("methods", ["auto"])[0],
        threads = lambda wildcards: config["resources"].get("tree", {}).get("threads", 4)
    threads: lambda wildcards: config["resources"].get("tree", {}).get("threads", 4)
    log:
        "logs/advanced_tree_{pathogen}.log"
    benchmark:
        "benchmarks/advanced_tree_{pathogen}.txt"
    resources:
        mem_mb = lambda wildcards: config["resources"].get("tree", {}).get("mem_mb", 8000),
        runtime = lambda wildcards: config["resources"].get("tree", {}).get("runtime", 120)
    shell:
        """
        python scripts/advanced_phylogenetics.py \
            {input.alignment} \
            {output.tree} \
            --method {params.method} \
            --threads {params.threads} \
            --metadata {input.metadata} \
            --stats-output {output.stats} \
            > {log} 2>&1
        """

# Temporal analysis and molecular clock estimation
rule temporal_analysis:
    input:
        tree = "results/tree/{pathogen}.nwk",
        metadata = get_input_metadata
    output:
        temporal_dir = directory("results/temporal/{pathogen}"),
        summary = "results/temporal/{pathogen}/temporal_analysis_summary.txt",
        results = "results/temporal/{pathogen}/temporal_analysis_results.json"
    params:
        threads = lambda wildcards: config["resources"].get("temporal", {}).get("threads", 2)
    threads: lambda wildcards: config["resources"].get("temporal", {}).get("threads", 2)
    log:
        "logs/temporal_analysis_{pathogen}.log"
    benchmark:
        "benchmarks/temporal_analysis_{pathogen}.txt"
    resources:
        mem_mb = lambda wildcards: config["resources"].get("temporal", {}).get("mem_mb", 4000)
    shell:
        """
        mkdir -p {output.temporal_dir}
        
        python scripts/temporal_analysis.py \
            {input.metadata} \
            --tree {input.tree} \
            --output-dir {output.temporal_dir} \
            > {log} 2>&1
        """

# Quality assessment of phylogenetic analysis
rule phylogenetic_qc:
    input:
        tree = "results/tree/{pathogen}.nwk",
        alignment = "results/aligned/{pathogen}_aligned.fasta",
        tree_stats = "results/tree/{pathogen}_tree_stats.json",
        temporal_results = "results/temporal/{pathogen}/temporal_analysis_results.json"
    output:
        qc_report = "results/qc_reports/{pathogen}_phylogenetic_qc.json",
        qc_summary = "results/qc_reports/{pathogen}_phylogenetic_qc.txt"
    log:
        "logs/phylogenetic_qc_{pathogen}.log"
    benchmark:
        "benchmarks/phylogenetic_qc_{pathogen}.txt"
    run:
        import json
        from datetime import datetime
        
        # Initialize QC results
        qc_results = {
            "timestamp": datetime.now().isoformat(),
            "pathogen": wildcards.pathogen,
            "tree_qc": {},
            "alignment_qc": {},
            "temporal_qc": {},
            "overall_quality": "unknown"
        }
        
        try:
            # Tree QC - basic checks
            import os
            qc_results["tree_qc"] = {
                "tree_exists": os.path.exists(input.tree),
                "tree_size": os.path.getsize(input.tree) if os.path.exists(input.tree) else 0
            }
            
            # Load tree statistics if available
            if os.path.exists(input.tree_stats):
                with open(input.tree_stats, 'r') as f:
                    tree_stats = json.load(f)
                qc_results["tree_qc"].update(tree_stats.get("tree_stats", {}))
            
            # Alignment QC
            qc_results["alignment_qc"] = {
                "alignment_exists": os.path.exists(input.alignment),
                "alignment_size": os.path.getsize(input.alignment) if os.path.exists(input.alignment) else 0
            }
            
            # Temporal QC
            if os.path.exists(input.temporal_results):
                with open(input.temporal_results, 'r') as f:
                    temporal_data = json.load(f)
                qc_results["temporal_qc"] = temporal_data.get("basic_temporal_stats", {})
            
            # Overall quality assessment
            issues = []
            
            if not qc_results["tree_qc"]["tree_exists"]:
                issues.append("Tree file not found")
                
            if not qc_results["alignment_qc"]["alignment_exists"]:
                issues.append("Alignment file not found")
                
            if qc_results["tree_qc"].get("num_terminals", 0) < 3:
                issues.append("Too few sequences for robust phylogenetic analysis")
                
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
            f.write("PHYLOGENETIC ANALYSIS QUALITY CONTROL REPORT\\n")
            f.write("=" * 50 + "\\n\\n")
            f.write(f"Pathogen: {wildcards.pathogen}\\n")
            f.write(f"Overall Quality: {qc_results['overall_quality'].upper()}\\n\\n")
            
            f.write("File Status:\\n")
            f.write(f"  - Tree: {'✓' if qc_results['tree_qc']['tree_exists'] else '✗'}\\n")
            f.write(f"  - Alignment: {'✓' if qc_results['alignment_qc']['alignment_exists'] else '✗'}\\n\\n")
                
            if qc_results.get("issues"):
                f.write("Issues Identified:\\n")
                for issue in qc_results["issues"]:
                    f.write(f"  - {issue}\\n")
                f.write("\\n")
                
            f.write(f"Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")

# Rule to generate advanced phylogenetic reports
rule advanced_phylogenetic_report:
    input:
        qc_report = "results/qc_reports/{pathogen}_phylogenetic_qc.json"
    output:
        html_report = "results/advanced_phylogenetics/{pathogen}_advanced_report.html",
        temporal_html = "results/advanced_phylogenetics/{pathogen}_temporal_analysis.html"
    log:
        "logs/advanced_phylogenetic_report_{pathogen}.log"
    run:
        import json
        from datetime import datetime
        
        # Load QC data
        with open(input.qc_report, 'r') as f:
            qc_data = json.load(f)
        
        # Generate HTML report
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Advanced Phylogenetic Analysis Report - {wildcards.pathogen}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                .status-good {{ color: green; font-weight: bold; }}
                .status-poor {{ color: red; font-weight: bold; }}
                .status-acceptable {{ color: orange; font-weight: bold; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Advanced Phylogenetic Analysis Report</h1>
                <p><strong>Pathogen:</strong> {wildcards.pathogen}</p>
                <p><strong>Analysis Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>Overall Quality:</strong> 
                    <span class="status-{qc_data.get('overall_quality', 'unknown')}">{qc_data.get('overall_quality', 'unknown').upper()}</span>
                </p>
            </div>
            
            <div class="section">
                <h2>File Status</h2>
                <table>
                    <tr><th>Component</th><th>Status</th></tr>
                    <tr><td>Phylogenetic Tree</td><td>{'✓ Present' if qc_data.get('tree_qc', {}).get('tree_exists') else '✗ Missing'}</td></tr>
                    <tr><td>Sequence Alignment</td><td>{'✓ Present' if qc_data.get('alignment_qc', {}).get('alignment_exists') else '✗ Missing'}</td></tr>
                </table>
            </div>
            
            <div class="section">
                <h2>Analysis Summary</h2>
                <p>This report summarizes the advanced phylogenetic analysis results for {wildcards.pathogen}.</p>
                <p>The analysis includes sequence alignment, tree building, and quality control assessment.</p>
            </div>
            
            <div class="section">
                <h2>Generated Files</h2>
                <ul>
                    <li>Tree file: results/tree/{wildcards.pathogen}.nwk</li>
                    <li>Alignment: results/aligned/{wildcards.pathogen}_aligned.fasta</li>
                    <li>Quality report: results/qc_reports/{wildcards.pathogen}_phylogenetic_qc.json</li>
                </ul>
            </div>
        </body>
        </html>
        """
        
        with open(output.html_report, 'w') as f:
            f.write(html_content)
        
        # Generate temporal analysis HTML
        temporal_html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Temporal Analysis Report - {wildcards.pathogen}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Temporal Analysis Report</h1>
                <p><strong>Pathogen:</strong> {wildcards.pathogen}</p>
                <p><strong>Analysis Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="section">
                <h2>Temporal Patterns</h2>
                <p>Temporal analysis results for {wildcards.pathogen} phylogenetic data.</p>
                <p>Results directory: results/temporal/{wildcards.pathogen}/</p>
            </div>
        </body>
        </html>
        """
        
        with open(output.temporal_html, 'w') as f:
            f.write(temporal_html)
