#!/usr/bin/env python3
"""
Enhanced QC Integration Script for RVF Nextstrain Pipeline

This script integrates the comprehensive QC system with the existing pipeline,
providing automated quality assessment and reporting capabilities.
"""

import os
import sys
import json
import argparse
import logging
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any

# Add scripts directory to path
sys.path.insert(0, "scripts")

try:
    from comprehensive_qc import ComprehensiveQC
except ImportError:
    logger = logging.getLogger(__name__)
    logger.error("Could not import comprehensive_qc module. Ensure it's in the scripts directory.")
    sys.exit(1)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class QCIntegration:
    """Quality Control Integration for RVF Nextstrain Pipeline"""
    
    def __init__(self, config_file: str = "config/qc_config.json"):
        """Initialize QC integration with configuration"""
        self.config_file = config_file
        self.config = self._load_config()
        self.results = {}
        
    def _load_config(self) -> Dict[str, Any]:
        """Load QC configuration"""
        default_config = {
            "quality_control": {
                "sequence_qc": {"min_length": 200, "max_n_content": 5.0},
                "metadata_qc": {"min_completeness_rate": 0.8},
                "nextclade_qc": {"min_qc_score": 80}
            },
            "thresholds": {
                "quality_scores": {"excellent": 85, "good": 70, "acceptable": 50}
            }
        }
        
        if os.path.exists(self.config_file):
            try:
                with open(self.config_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Could not load config from {self.config_file}: {e}")
                logger.warning("Using default configuration")
        
        return default_config
    
    def run_pipeline_qc(self, stage: str = "all") -> Dict[str, Any]:
        """
        Run quality control for different pipeline stages
        
        Args:
            stage: QC stage to run ("download", "filter", "align", "tree", "all")
            
        Returns:
            Dictionary containing QC results
        """
        logger.info(f"Running QC for stage: {stage}")
        
        if stage == "all" or stage == "download":
            self.results["download"] = self._qc_download_stage()
        
        if stage == "all" or stage == "filter":
            self.results["filter"] = self._qc_filter_stage()
        
        if stage == "all" or stage == "align":
            self.results["align"] = self._qc_align_stage()
        
        if stage == "all" or stage == "tree":
            self.results["tree"] = self._qc_tree_stage()
        
        # Generate comprehensive report
        self._generate_integration_report()
        
        return self.results
    
    def _qc_download_stage(self) -> Dict[str, Any]:
        """Quality control for download stage"""
        logger.info("Running QC for download stage...")
        
        qc_results = {
            "stage": "download",
            "timestamp": datetime.now().isoformat(),
            "status": "unknown",
            "metrics": {},
            "issues": [],
            "recommendations": []
        }
        
        # Check downloaded data files
        raw_metadata = "data/metadata/raw/rvf_metadata.tsv"
        raw_sequences = "data/sequences/raw/rvf_sequences.fasta"
        
        if os.path.exists(raw_metadata) and os.path.exists(raw_sequences):
            # Run comprehensive QC on raw data
            qc_system = ComprehensiveQC(self.config)
            
            try:
                results = qc_system.run_comprehensive_qc(
                    metadata_file=raw_metadata,
                    sequences_file=raw_sequences,
                    output_dir="results/qc_reports/download"
                )
                
                qc_results["metrics"] = {
                    "metadata_score": results.get("metadata_qc", {}).get("quality_score", 0),
                    "sequence_score": results.get("sequence_qc", {}).get("quality_score", 0),
                    "overall_score": results.get("summary", {}).get("overall_score", 0)
                }
                
                qc_results["status"] = results.get("overall_quality", "unknown")
                qc_results["recommendations"] = results.get("recommendations", [])
                
            except Exception as e:
                logger.error(f"Error in download stage QC: {e}")
                qc_results["status"] = "error"
                qc_results["issues"].append(str(e))
        else:
            qc_results["status"] = "missing_files"
            qc_results["issues"].append("Raw data files not found")
        
        return qc_results
    
    def _qc_filter_stage(self) -> Dict[str, Any]:
        """Quality control for filter stage"""
        logger.info("Running QC for filter stage...")
        
        qc_results = {
            "stage": "filter",
            "timestamp": datetime.now().isoformat(),
            "status": "unknown",
            "metrics": {},
            "filtering_stats": {},
            "issues": [],
            "recommendations": []
        }
        
        # Check filtered data files
        filtered_metadata = "results/filtered/rvf_metadata.tsv"
        filtered_sequences = "results/filtered/rvf_filtered.fasta"
        
        if os.path.exists(filtered_metadata) and os.path.exists(filtered_sequences):
            # Run comprehensive QC on filtered data
            qc_system = ComprehensiveQC(self.config)
            
            try:
                results = qc_system.run_comprehensive_qc(
                    metadata_file=filtered_metadata,
                    sequences_file=filtered_sequences,
                    output_dir="results/qc_reports/filter"
                )
                
                # Calculate filtering statistics
                if os.path.exists("data/metadata/raw/rvf_metadata.tsv"):
                    import pandas as pd
                    raw_metadata = pd.read_csv("data/metadata/raw/rvf_metadata.tsv", sep='\t')
                    filtered_metadata_df = pd.read_csv(filtered_metadata, sep='\t')
                    
                    qc_results["filtering_stats"] = {
                        "input_sequences": len(raw_metadata),
                        "output_sequences": len(filtered_metadata_df),
                        "retention_rate": len(filtered_metadata_df) / len(raw_metadata),
                        "filtered_out": len(raw_metadata) - len(filtered_metadata_df)
                    }
                
                qc_results["metrics"] = {
                    "metadata_score": results.get("metadata_qc", {}).get("quality_score", 0),
                    "sequence_score": results.get("sequence_qc", {}).get("quality_score", 0),
                    "geographic_score": results.get("geographic_qc", {}).get("quality_score", 0),
                    "overall_score": results.get("summary", {}).get("overall_score", 0)
                }
                
                qc_results["status"] = results.get("overall_quality", "unknown")
                qc_results["recommendations"] = results.get("recommendations", [])
                
                # Check if retention rate is acceptable
                retention_rate = qc_results["filtering_stats"].get("retention_rate", 0)
                min_retention = self.config.get("thresholds", {}).get("data_retention", {}).get("min_retention_rate", 0.3)
                
                if retention_rate < min_retention:
                    qc_results["issues"].append(f"Low retention rate: {retention_rate:.2%} < {min_retention:.2%}")
                    qc_results["recommendations"].append("Review filtering criteria - retention rate is too low")
                
            except Exception as e:
                logger.error(f"Error in filter stage QC: {e}")
                qc_results["status"] = "error"
                qc_results["issues"].append(str(e))
        else:
            qc_results["status"] = "missing_files"
            qc_results["issues"].append("Filtered data files not found")
        
        return qc_results
    
    def _qc_align_stage(self) -> Dict[str, Any]:
        """Quality control for alignment stage"""
        logger.info("Running QC for alignment stage...")
        
        qc_results = {
            "stage": "align",
            "timestamp": datetime.now().isoformat(),
            "status": "unknown",
            "metrics": {},
            "alignment_stats": {},
            "issues": [],
            "recommendations": []
        }
        
        # Check alignment files
        alignment_file = "results/aligned/rvf_aligned.fasta"
        
        if os.path.exists(alignment_file):
            try:
                from Bio import AlignIO
                
                # Load alignment
                alignment = AlignIO.read(alignment_file, "fasta")
                
                qc_results["alignment_stats"] = {
                    "sequence_count": len(alignment),
                    "alignment_length": alignment.get_alignment_length(),
                    "conserved_positions": self._count_conserved_positions(alignment),
                    "gap_statistics": self._calculate_gap_statistics(alignment)
                }
                
                # Assess alignment quality
                alignment_length = qc_results["alignment_stats"]["alignment_length"]
                gap_rate = qc_results["alignment_stats"]["gap_statistics"]["mean_gap_rate"]
                
                if gap_rate > 0.5:
                    qc_results["issues"].append(f"High gap rate in alignment: {gap_rate:.2%}")
                    qc_results["status"] = "poor"
                elif gap_rate > 0.3:
                    qc_results["status"] = "acceptable"
                else:
                    qc_results["status"] = "good"
                
                qc_results["metrics"]["alignment_score"] = max(0, 100 * (1 - gap_rate))
                
            except Exception as e:
                logger.error(f"Error in alignment stage QC: {e}")
                qc_results["status"] = "error"
                qc_results["issues"].append(str(e))
        else:
            qc_results["status"] = "missing_files"
            qc_results["issues"].append("Alignment file not found")
        
        return qc_results
    
    def _qc_tree_stage(self) -> Dict[str, Any]:
        """Quality control for tree stage"""
        logger.info("Running QC for tree stage...")
        
        qc_results = {
            "stage": "tree",
            "timestamp": datetime.now().isoformat(),
            "status": "unknown",
            "metrics": {},
            "tree_stats": {},
            "issues": [],
            "recommendations": []
        }
        
        # Check tree files
        tree_file = "results/tree/rvf.nwk"
        metadata_file = "results/filtered/rvf_metadata.tsv"
        
        if os.path.exists(tree_file):
            try:
                # Run phylogenetic QC
                qc_system = ComprehensiveQC(self.config)
                results = qc_system._phylogenetic_qc(tree_file, 
                    pd.read_csv(metadata_file, sep='\t') if os.path.exists(metadata_file) else pd.DataFrame())
                
                qc_results["tree_stats"] = results.get("tree_statistics", {})
                qc_results["metrics"]["phylogenetic_score"] = results.get("quality_score", 0)
                
                # Assess tree quality
                terminal_nodes = qc_results["tree_stats"].get("terminal_nodes", 0)
                min_sequences = self.config.get("quality_control", {}).get("phylogenetic_qc", {}).get("min_sequences", 10)
                
                if terminal_nodes < min_sequences:
                    qc_results["issues"].append(f"Too few sequences in tree: {terminal_nodes} < {min_sequences}")
                    qc_results["status"] = "poor"
                elif terminal_nodes < min_sequences * 2:
                    qc_results["status"] = "acceptable"
                else:
                    qc_results["status"] = "good"
                
            except Exception as e:
                logger.error(f"Error in tree stage QC: {e}")
                qc_results["status"] = "error"
                qc_results["issues"].append(str(e))
        else:
            qc_results["status"] = "missing_files"
            qc_results["issues"].append("Tree file not found")
        
        return qc_results
    
    def _count_conserved_positions(self, alignment) -> int:
        """Count conserved positions in alignment"""
        conserved = 0
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            if len(set(column.upper().replace('-', '').replace('N', ''))) == 1:
                conserved += 1
        return conserved
    
    def _calculate_gap_statistics(self, alignment) -> Dict[str, float]:
        """Calculate gap statistics for alignment"""
        gap_rates = []
        for record in alignment:
            seq_str = str(record.seq)
            gap_rate = seq_str.count('-') / len(seq_str)
            gap_rates.append(gap_rate)
        
        return {
            "mean_gap_rate": np.mean(gap_rates) if gap_rates else 0,
            "max_gap_rate": max(gap_rates) if gap_rates else 0,
            "min_gap_rate": min(gap_rates) if gap_rates else 0
        }
    
    def _generate_integration_report(self):
        """Generate comprehensive integration report"""
        logger.info("Generating QC integration report...")
        
        # Create QC reports directory
        report_dir = "results/qc_reports/integration"
        os.makedirs(report_dir, exist_ok=True)
        
        # Generate JSON report
        report_file = os.path.join(report_dir, "pipeline_qc_integration.json")
        with open(report_file, 'w') as f:
            json.dump(self.results, f, indent=2, default=str)
        
        # Generate summary report
        summary_file = os.path.join(report_dir, "pipeline_qc_summary.txt")
        self._generate_summary_report(summary_file)
        
        # Generate HTML dashboard
        html_file = os.path.join(report_dir, "pipeline_qc_dashboard.html")
        self._generate_html_dashboard(html_file)
        
        logger.info(f"QC integration reports generated in {report_dir}")
    
    def _generate_summary_report(self, output_file: str):
        """Generate text summary report"""
        with open(output_file, 'w') as f:
            f.write("RVF NEXTSTRAIN PIPELINE - QUALITY CONTROL INTEGRATION SUMMARY\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Report Generated: {datetime.now().isoformat()}\n\n")
            
            # Overall pipeline status
            statuses = [result.get("status", "unknown") for result in self.results.values()]
            if all(s in ["good", "excellent"] for s in statuses):
                overall_status = "EXCELLENT"
            elif all(s in ["good", "excellent", "acceptable"] for s in statuses):
                overall_status = "GOOD"
            elif any(s == "error" for s in statuses):
                overall_status = "ERROR"
            else:
                overall_status = "NEEDS ATTENTION"
            
            f.write(f"Overall Pipeline Status: {overall_status}\n\n")
            
            # Stage-by-stage results
            f.write("STAGE-BY-STAGE RESULTS:\n")
            f.write("-" * 30 + "\n")
            
            for stage, result in self.results.items():
                f.write(f"\n{stage.upper()} STAGE:\n")
                f.write(f"  Status: {result.get('status', 'unknown').upper()}\n")
                
                if "metrics" in result:
                    f.write("  Metrics:\n")
                    for metric, value in result["metrics"].items():
                        f.write(f"    {metric}: {value}\n")
                
                if result.get("issues"):
                    f.write("  Issues:\n")
                    for issue in result["issues"]:
                        f.write(f"    - {issue}\n")
                
                if result.get("recommendations"):
                    f.write("  Recommendations:\n")
                    for rec in result["recommendations"]:
                        f.write(f"    - {rec}\n")
    
    def _generate_html_dashboard(self, output_file: str):
        """Generate HTML dashboard"""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>RVF Nextstrain Pipeline - QC Dashboard</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f8ff; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
        .stage-card {{ border: 1px solid #ddd; border-radius: 5px; margin: 10px 0; padding: 15px; }}
        .status-excellent {{ border-left: 5px solid #28a745; }}
        .status-good {{ border-left: 5px solid #007bff; }}
        .status-acceptable {{ border-left: 5px solid #ffc107; }}
        .status-poor {{ border-left: 5px solid #dc3545; }}
        .status-error {{ border-left: 5px solid #dc3545; background-color: #fff5f5; }}
        .metrics {{ display: flex; flex-wrap: wrap; }}
        .metric {{ margin: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 3px; }}
        .issues {{ color: #dc3545; }}
        .recommendations {{ color: #007bff; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>RVF Nextstrain Pipeline - Quality Control Dashboard</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    {self._generate_stage_cards_html()}
    
    <div class="footer">
        <p><em>This dashboard provides real-time quality control monitoring for the RVF Nextstrain pipeline.</em></p>
    </div>
</body>
</html>
        """
        
        with open(output_file, 'w') as f:
            f.write(html_content)
    
    def _generate_stage_cards_html(self) -> str:
        """Generate HTML cards for each stage"""
        cards_html = ""
        
        for stage, result in self.results.items():
            status = result.get("status", "unknown")
            
            cards_html += f"""
            <div class="stage-card status-{status}">
                <h3>{stage.upper()} Stage</h3>
                <p><strong>Status:</strong> {status.upper()}</p>
                
                {self._generate_metrics_html(result.get("metrics", {}))}
                
                {self._generate_issues_html(result.get("issues", []))}
                
                {self._generate_recommendations_html(result.get("recommendations", []))}
            </div>
            """
        
        return cards_html
    
    def _generate_metrics_html(self, metrics: Dict[str, Any]) -> str:
        """Generate metrics HTML"""
        if not metrics:
            return ""
        
        html = "<div class='metrics'>"
        for metric, value in metrics.items():
            if isinstance(value, float):
                html += f"<div class='metric'><strong>{metric}:</strong> {value:.1f}</div>"
            else:
                html += f"<div class='metric'><strong>{metric}:</strong> {value}</div>"
        html += "</div>"
        
        return html
    
    def _generate_issues_html(self, issues: List[str]) -> str:
        """Generate issues HTML"""
        if not issues:
            return ""
        
        html = "<div class='issues'><strong>Issues:</strong><ul>"
        for issue in issues:
            html += f"<li>{issue}</li>"
        html += "</ul></div>"
        
        return html
    
    def _generate_recommendations_html(self, recommendations: List[str]) -> str:
        """Generate recommendations HTML"""
        if not recommendations:
            return ""
        
        html = "<div class='recommendations'><strong>Recommendations:</strong><ul>"
        for rec in recommendations:
            html += f"<li>{rec}</li>"
        html += "</ul></div>"
        
        return html
    
    def run_nextclade_qc(self, segment: str = "all") -> Dict[str, Any]:
        """Run Nextclade quality control"""
        logger.info(f"Running Nextclade QC for segment: {segment}")
        
        # Check if Nextclade files exist
        nextclade_dir = "results/nextclade"
        if not os.path.exists(nextclade_dir):
            logger.warning("Nextclade results directory not found")
            return {"status": "not_run", "message": "Nextclade not executed"}
        
        # Run Nextclade if not already done
        self._ensure_nextclade_run()
        
        # Analyze Nextclade results
        nextclade_json = os.path.join(nextclade_dir, "rvf_nextclade.json")
        if os.path.exists(nextclade_json):
            qc_system = ComprehensiveQC(self.config)
            return qc_system._nextclade_qc(nextclade_json)
        else:
            return {"status": "missing_results", "message": "Nextclade results not found"}
    
    def _ensure_nextclade_run(self):
        """Ensure Nextclade has been run"""
        nextclade_script = "scripts/run_nextclade_windows.py"
        if os.path.exists(nextclade_script):
            try:
                subprocess.run([sys.executable, nextclade_script], check=True)
                logger.info("Nextclade QC completed")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running Nextclade: {e}")
        else:
            logger.warning("Nextclade script not found")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="QC Integration for RVF Nextstrain Pipeline")
    parser.add_argument("--stage", choices=["download", "filter", "align", "tree", "all"], 
                       default="all", help="Pipeline stage to run QC for")
    parser.add_argument("--config", default="config/qc_config.json", 
                       help="QC configuration file")
    parser.add_argument("--output-dir", default="results/qc_reports", 
                       help="Output directory for reports")
    parser.add_argument("--run-nextclade", action="store_true", 
                       help="Run Nextclade QC")
    
    args = parser.parse_args()
    
    # Initialize QC integration
    qc_integration = QCIntegration(args.config)
    
    # Run pipeline QC
    results = qc_integration.run_pipeline_qc(args.stage)
    
    # Run Nextclade QC if requested
    if args.run_nextclade:
        nextclade_results = qc_integration.run_nextclade_qc()
        results["nextclade"] = nextclade_results
    
    # Print summary
    print(f"\nQuality Control Integration Complete!")
    print(f"Stage: {args.stage}")
    
    for stage, result in results.items():
        status = result.get("status", "unknown")
        print(f"{stage}: {status.upper()}")
    
    print(f"Reports available in: {args.output_dir}")


if __name__ == "__main__":
    # Add numpy import for gap statistics calculation
    try:
        import numpy as np
        import pandas as pd
    except ImportError:
        print("Required packages not found. Please install numpy and pandas.")
        sys.exit(1)
    
    main()
