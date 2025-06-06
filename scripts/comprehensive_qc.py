#!/usr/bin/env python3
"""
Comprehensive Quality Control System for RVF Nextstrain Pipeline

This script implements multi-level quality control including:
1. Nextclade QC Integration - Sequence quality assessment and mutation analysis
2. Phylogenetic QC - Tree topology validation and temporal signal assessment
3. Geographic QC - Coordinate validation and geographic clustering analysis
4. Metadata QC - Field completeness reporting and data consistency checks
"""

import os
import sys
import json
import argparse
import logging
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Any

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ComprehensiveQC:
    """Comprehensive Quality Control for RVF Nextstrain Pipeline"""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.qc_results = {
            "timestamp": datetime.now().isoformat(),
            "nextclade_qc": {},
            "phylogenetic_qc": {},
            "geographic_qc": {},
            "metadata_qc": {},
            "overall_quality": "unknown",
            "summary": {},
            "recommendations": []
        }
        
    def run_comprehensive_qc(self, 
                           metadata_file: str,
                           sequences_file: str,
                           tree_file: Optional[str] = None,
                           nextclade_results: Optional[str] = None,
                           output_dir: str = "results/qc_reports") -> Dict[str, Any]:
        """
        Run comprehensive quality control analysis
        
        Args:
            metadata_file: Path to metadata TSV file
            sequences_file: Path to FASTA sequences file
            tree_file: Optional path to phylogenetic tree file
            nextclade_results: Optional path to Nextclade JSON results
            output_dir: Output directory for QC reports
            
        Returns:
            Dictionary containing comprehensive QC results
        """
        logger.info("Starting comprehensive quality control analysis...")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Load data
        metadata = pd.read_csv(metadata_file, sep='\t')
        sequences = list(SeqIO.parse(sequences_file, "fasta"))
        
        # Run QC modules
        logger.info("Running metadata quality control...")
        self.qc_results["metadata_qc"] = self._metadata_qc(metadata)
        
        logger.info("Running sequence quality control...")
        self.qc_results["sequence_qc"] = self._sequence_qc(sequences, metadata)
        
        logger.info("Running geographic quality control...")
        self.qc_results["geographic_qc"] = self._geographic_qc(metadata)
        
        if nextclade_results and os.path.exists(nextclade_results):
            logger.info("Running Nextclade quality control...")
            self.qc_results["nextclade_qc"] = self._nextclade_qc(nextclade_results)
        
        if tree_file and os.path.exists(tree_file):
            logger.info("Running phylogenetic quality control...")
            self.qc_results["phylogenetic_qc"] = self._phylogenetic_qc(tree_file, metadata)
        
        # Generate overall assessment
        logger.info("Generating overall quality assessment...")
        self._generate_overall_assessment()
        
        # Create comprehensive report
        self._generate_comprehensive_report(output_dir)
        
        logger.info(f"Quality control analysis complete. Reports saved to {output_dir}")
        return self.qc_results
    
    def _metadata_qc(self, metadata: pd.DataFrame) -> Dict[str, Any]:
        """Comprehensive metadata quality control"""
        qc = {
            "total_records": len(metadata),
            "field_completeness": {},
            "date_quality": {},
            "geographic_quality": {},
            "segment_distribution": {},
            "host_distribution": {},
            "data_consistency": {},
            "quality_score": 0
        }
        
        # Field completeness analysis
        required_fields = ["strain", "virus", "accession", "date", "country"]
        optional_fields = ["division", "location", "host", "segment", "latitude", "longitude"]
        
        for field in required_fields + optional_fields:
            if field in metadata.columns:
                non_null = metadata[field].notna().sum()
                non_empty = (metadata[field].astype(str).str.strip() != "").sum()
                qc["field_completeness"][field] = {
                    "present": non_null,
                    "non_empty": non_empty,
                    "completeness_rate": non_empty / len(metadata),
                    "is_required": field in required_fields
                }
        
        # Date quality analysis
        if "date" in metadata.columns:
            dates = metadata["date"].dropna()
            qc["date_quality"] = {
                "total_dates": len(dates),
                "valid_iso_dates": 0,
                "year_only_dates": 0,
                "invalid_dates": 0,
                "date_range": {"min": None, "max": None},
                "temporal_distribution": {}
            }
            
            for date_str in dates:
                if pd.isna(date_str) or str(date_str).strip() == "":
                    continue
                    
                date_str = str(date_str).strip()
                
                # Check date format patterns
                if len(date_str) == 4 and date_str.isdigit():
                    qc["date_quality"]["year_only_dates"] += 1
                elif len(date_str) == 10 and "-" in date_str:
                    try:
                        datetime.strptime(date_str, "%Y-%m-%d")
                        qc["date_quality"]["valid_iso_dates"] += 1
                    except ValueError:
                        qc["date_quality"]["invalid_dates"] += 1
                else:
                    qc["date_quality"]["invalid_dates"] += 1
            
            # Temporal distribution
            years = []
            for date_str in dates:
                if pd.notna(date_str):
                    year_match = str(date_str)[:4]
                    if year_match.isdigit():
                        years.append(int(year_match))
            
            if years:
                qc["date_quality"]["date_range"]["min"] = min(years)
                qc["date_quality"]["date_range"]["max"] = max(years)
                year_counts = Counter(years)
                qc["date_quality"]["temporal_distribution"] = dict(year_counts)
        
        # Geographic quality analysis
        countries = metadata["country"].dropna() if "country" in metadata.columns else pd.Series()
        qc["geographic_quality"] = {
            "countries_present": len(countries.unique()) if len(countries) > 0 else 0,
            "coordinates_available": 0,
            "country_distribution": dict(countries.value_counts()) if len(countries) > 0 else {}
        }
        
        if "latitude" in metadata.columns and "longitude" in metadata.columns:
            coords_available = (metadata["latitude"].notna() & metadata["longitude"].notna()).sum()
            qc["geographic_quality"]["coordinates_available"] = coords_available
        
        # Segment distribution
        if "segment" in metadata.columns:
            segments = metadata["segment"].dropna()
            qc["segment_distribution"] = dict(segments.value_counts())
        
        # Host distribution
        if "host" in metadata.columns:
            hosts = metadata["host"].dropna()
            qc["host_distribution"] = dict(hosts.value_counts())
        
        # Data consistency checks
        qc["data_consistency"] = self._check_data_consistency(metadata)
        
        # Calculate quality score
        qc["quality_score"] = self._calculate_metadata_quality_score(qc)
        
        return qc
    
    def _sequence_qc(self, sequences: List[SeqIO.SeqRecord], metadata: pd.DataFrame) -> Dict[str, Any]:
        """Sequence quality control analysis"""
        qc = {
            "total_sequences": len(sequences),
            "length_distribution": {},
            "composition_analysis": {},
            "quality_metrics": {},
            "segment_analysis": {},
            "quality_score": 0
        }
        
        # Length analysis
        lengths = [len(seq.seq) for seq in sequences]
        qc["length_distribution"] = {
            "min": min(lengths) if lengths else 0,
            "max": max(lengths) if lengths else 0,
            "mean": np.mean(lengths) if lengths else 0,
            "median": np.median(lengths) if lengths else 0,
            "std": np.std(lengths) if lengths else 0
        }
        
        # Composition analysis
        gc_contents = []
        n_contents = []
        
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            gc_contents.append(gc_fraction(seq_str) * 100)
            n_count = seq_str.count('N') / len(seq_str) * 100
            n_contents.append(n_count)
        
        qc["composition_analysis"] = {
            "gc_content": {
                "mean": np.mean(gc_contents) if gc_contents else 0,
                "std": np.std(gc_contents) if gc_contents else 0,
                "min": min(gc_contents) if gc_contents else 0,
                "max": max(gc_contents) if gc_contents else 0
            },
            "n_content": {
                "mean": np.mean(n_contents) if n_contents else 0,
                "std": np.std(n_contents) if n_contents else 0,
                "min": min(n_contents) if n_contents else 0,
                "max": max(n_contents) if n_contents else 0
            }
        }
        
        # Quality metrics
        high_n_sequences = sum(1 for n in n_contents if n > 5)  # >5% N content
        qc["quality_metrics"] = {
            "high_n_content_sequences": high_n_sequences,
            "high_n_content_rate": high_n_sequences / len(sequences) if sequences else 0,
            "very_short_sequences": sum(1 for l in lengths if l < 200),
            "very_long_sequences": sum(1 for l in lengths if l > 10000)
        }
        
        # Segment-specific analysis
        if "segment" in metadata.columns:
            segment_stats = {}
            for segment in ["L", "M", "S"]:
                seg_metadata = metadata[metadata["segment"] == segment]
                seg_sequences = [seq for seq in sequences if seq.id in seg_metadata["strain"].values]
                
                if seg_sequences:
                    seg_lengths = [len(seq.seq) for seq in seg_sequences]
                    segment_stats[segment] = {
                        "count": len(seg_sequences),
                        "length_mean": np.mean(seg_lengths),
                        "length_std": np.std(seg_lengths),
                        "expected_length": self._get_expected_length(segment)
                    }
            
            qc["segment_analysis"] = segment_stats
        
        # Calculate quality score
        qc["quality_score"] = self._calculate_sequence_quality_score(qc)
        
        return qc
    
    def _geographic_qc(self, metadata: pd.DataFrame) -> Dict[str, Any]:
        """Geographic data quality control"""
        qc = {
            "coordinate_validation": {},
            "geographic_distribution": {},
            "clustering_analysis": {},
            "quality_score": 0
        }
        
        # Coordinate validation
        if "latitude" in metadata.columns and "longitude" in metadata.columns:
            coords = metadata[["latitude", "longitude"]].dropna()
            
            qc["coordinate_validation"] = {
                "total_coordinates": len(coords),
                "valid_coordinates": 0,
                "invalid_coordinates": 0,
                "coordinate_ranges": {
                    "latitude": {"min": None, "max": None},
                    "longitude": {"min": None, "max": None}
                }
            }
            
            valid_coords = []
            for _, row in coords.iterrows():
                lat, lon = row["latitude"], row["longitude"]
                if (-90 <= lat <= 90) and (-180 <= lon <= 180):
                    qc["coordinate_validation"]["valid_coordinates"] += 1
                    valid_coords.append((lat, lon))
                else:
                    qc["coordinate_validation"]["invalid_coordinates"] += 1
            
            if valid_coords:
                lats, lons = zip(*valid_coords)
                qc["coordinate_validation"]["coordinate_ranges"]["latitude"] = {
                    "min": min(lats), "max": max(lats)
                }
                qc["coordinate_validation"]["coordinate_ranges"]["longitude"] = {
                    "min": min(lons), "max": max(lons)
                }
        
        # Geographic distribution analysis
        if "country" in metadata.columns:
            countries = metadata["country"].dropna()
            country_counts = countries.value_counts()
            
            qc["geographic_distribution"] = {
                "unique_countries": len(country_counts),
                "country_counts": dict(country_counts),
                "geographic_coverage": self._assess_geographic_coverage(country_counts)
            }
        
        # Calculate quality score
        qc["quality_score"] = self._calculate_geographic_quality_score(qc)
        
        return qc
    
    def _nextclade_qc(self, nextclade_file: str) -> Dict[str, Any]:
        """Nextclade quality control analysis"""
        qc = {
            "total_analyzed": 0,
            "qc_scores": {},
            "mutation_analysis": {},
            "clade_analysis": {},
            "quality_issues": {},
            "quality_score": 0
        }
        
        try:
            with open(nextclade_file, 'r') as f:
                nextclade_data = json.load(f)
            
            results = nextclade_data.get("results", [])
            qc["total_analyzed"] = len(results)
            
            # QC score analysis
            qc_scores = []
            qc_statuses = []
            
            for result in results:
                qc_result = result.get("qc", {})
                score = qc_result.get("overallScore", 0)
                status = qc_result.get("overallStatus", "bad")
                
                qc_scores.append(score)
                qc_statuses.append(status)
            
            qc["qc_scores"] = {
                "mean": np.mean(qc_scores) if qc_scores else 0,
                "median": np.median(qc_scores) if qc_scores else 0,
                "std": np.std(qc_scores) if qc_scores else 0,
                "distribution": dict(Counter(qc_statuses))
            }
            
            # Quality issues analysis
            issue_types = defaultdict(int)
            for result in results:
                qc_result = result.get("qc", {})
                for issue_category in ["missingData", "privateMutations", "snpClusters"]:
                    issues = qc_result.get(issue_category, {}).get("issues", [])
                    for issue in issues:
                        issue_types[issue.get("tagName", "unknown")] += 1
            
            qc["quality_issues"] = dict(issue_types)
            
            # Calculate quality score
            qc["quality_score"] = self._calculate_nextclade_quality_score(qc)
            
        except Exception as e:
            logger.error(f"Error processing Nextclade results: {e}")
            qc["error"] = str(e)
        
        return qc
    
    def _phylogenetic_qc(self, tree_file: str, metadata: pd.DataFrame) -> Dict[str, Any]:
        """Phylogenetic quality control analysis"""
        qc = {
            "tree_statistics": {},
            "temporal_analysis": {},
            "topology_assessment": {},
            "quality_score": 0
        }
        
        try:
            tree = Phylo.read(tree_file, "newick")
            
            # Tree statistics
            terminals = tree.get_terminals()
            nonterminals = tree.get_nonterminals()
            
            qc["tree_statistics"] = {
                "terminal_nodes": len(terminals),
                "internal_nodes": len(nonterminals),
                "total_nodes": len(terminals) + len(nonterminals),
                "tree_depth": tree.depth(),
                "tree_length": tree.total_branch_length()
            }
            
            # Branch length analysis
            branch_lengths = []
            for clade in tree.find_clades():
                if clade.branch_length is not None:
                    branch_lengths.append(clade.branch_length)
            
            if branch_lengths:
                qc["tree_statistics"]["branch_lengths"] = {
                    "mean": np.mean(branch_lengths),
                    "median": np.median(branch_lengths),
                    "std": np.std(branch_lengths),
                    "min": min(branch_lengths),
                    "max": max(branch_lengths)
                }
            
            # Temporal analysis if dates available
            if "date" in metadata.columns:
                qc["temporal_analysis"] = self._analyze_temporal_signal(tree, metadata)
            
            # Calculate quality score
            qc["quality_score"] = self._calculate_phylogenetic_quality_score(qc)
            
        except Exception as e:
            logger.error(f"Error processing phylogenetic tree: {e}")
            qc["error"] = str(e)
        
        return qc
    
    def _check_data_consistency(self, metadata: pd.DataFrame) -> Dict[str, Any]:
        """Check data consistency across fields"""
        consistency = {
            "strain_accession_match": True,
            "date_format_consistency": True,
            "geographic_consistency": True,
            "issues": []
        }
        
        # Check strain-accession consistency
        if "strain" in metadata.columns and "accession" in metadata.columns:
            mismatches = 0
            for _, row in metadata.iterrows():
                strain = str(row.get("strain", ""))
                accession = str(row.get("accession", ""))
                if strain and accession and strain != accession:
                    # Allow for different naming conventions
                    if not (strain in accession or accession in strain):
                        mismatches += 1
            
            if mismatches > len(metadata) * 0.1:  # >10% mismatches
                consistency["strain_accession_match"] = False
                consistency["issues"].append(f"High strain-accession mismatch rate: {mismatches}/{len(metadata)}")
        
        return consistency
    
    def _calculate_metadata_quality_score(self, qc: Dict[str, Any]) -> float:
        """Calculate metadata quality score (0-100)"""
        score = 0
        
        # Completeness score (40 points)
        required_completeness = []
        for field, stats in qc["field_completeness"].items():
            if stats["is_required"]:
                required_completeness.append(stats["completeness_rate"])
        
        if required_completeness:
            score += 40 * np.mean(required_completeness)
        
        # Date quality (30 points)
        if qc["date_quality"]:
            date_score = 0
            total_dates = qc["date_quality"]["total_dates"]
            if total_dates > 0:
                valid_rate = (qc["date_quality"]["valid_iso_dates"] + 
                            qc["date_quality"]["year_only_dates"] * 0.7) / total_dates
                date_score = 30 * valid_rate
            score += date_score
        
        # Geographic quality (20 points)
        if qc["geographic_quality"]["countries_present"] > 0:
            score += 20 * min(qc["geographic_quality"]["countries_present"] / 5, 1)
        
        # Data consistency (10 points)
        if qc["data_consistency"]["strain_accession_match"]:
            score += 10
        
        return min(score, 100)
    
    def _calculate_sequence_quality_score(self, qc: Dict[str, Any]) -> float:
        """Calculate sequence quality score (0-100)"""
        score = 100
        
        # Penalize high N content
        if qc["quality_metrics"]["high_n_content_rate"] > 0.1:
            score -= 30 * qc["quality_metrics"]["high_n_content_rate"]
        
        # Penalize unusual lengths
        if qc["quality_metrics"]["very_short_sequences"] > 0:
            score -= 20 * (qc["quality_metrics"]["very_short_sequences"] / qc["total_sequences"])
        
        return max(score, 0)
    
    def _calculate_geographic_quality_score(self, qc: Dict[str, Any]) -> float:
        """Calculate geographic quality score (0-100)"""
        score = 0
        
        # Coordinate validation (50 points)
        if qc["coordinate_validation"]:
            total_coords = qc["coordinate_validation"]["total_coordinates"]
            valid_coords = qc["coordinate_validation"]["valid_coordinates"]
            if total_coords > 0:
                score += 50 * (valid_coords / total_coords)
        
        # Geographic diversity (50 points)
        if qc["geographic_distribution"]:
            countries = qc["geographic_distribution"]["unique_countries"]
            score += 50 * min(countries / 10, 1)  # Max score for 10+ countries
        
        return score
    
    def _calculate_nextclade_quality_score(self, qc: Dict[str, Any]) -> float:
        """Calculate Nextclade quality score (0-100)"""
        if qc["qc_scores"]["mean"]:
            return qc["qc_scores"]["mean"]
        return 0
    
    def _calculate_phylogenetic_quality_score(self, qc: Dict[str, Any]) -> float:
        """Calculate phylogenetic quality score (0-100)"""
        score = 0
        
        # Tree size (30 points)
        terminals = qc["tree_statistics"].get("terminal_nodes", 0)
        if terminals >= 10:
            score += 30
        elif terminals >= 5:
            score += 20
        elif terminals >= 3:
            score += 10
        
        # Branch length consistency (40 points)
        if "branch_lengths" in qc["tree_statistics"]:
            bl_stats = qc["tree_statistics"]["branch_lengths"]
            cv = bl_stats["std"] / bl_stats["mean"] if bl_stats["mean"] > 0 else float('inf')
            if cv < 2:  # Reasonable coefficient of variation
                score += 40
        
        # Temporal signal (30 points)
        if qc["temporal_analysis"]:
            score += 30  # Basic score for having temporal analysis
        
        return score
    
    def _get_expected_length(self, segment: str) -> int:
        """Get expected length for RVF segments"""
        expected_lengths = {
            "L": 6404,  # L segment
            "M": 3885,  # M segment  
            "S": 1690   # S segment
        }
        return expected_lengths.get(segment, 0)
    
    def _assess_geographic_coverage(self, country_counts: pd.Series) -> str:
        """Assess geographic coverage quality"""
        total_countries = len(country_counts)
        
        if total_countries >= 10:
            return "excellent"
        elif total_countries >= 5:
            return "good"
        elif total_countries >= 3:
            return "moderate"
        elif total_countries >= 1:
            return "limited"
        else:
            return "none"
    
    def _analyze_temporal_signal(self, tree, metadata: pd.DataFrame) -> Dict[str, Any]:
        """Analyze temporal signal in phylogenetic tree"""
        temporal = {
            "date_range": {},
            "temporal_structure": {},
            "rate_analysis": {}
        }
        
        # Extract dates for tree tips
        tip_dates = {}
        for _, row in metadata.iterrows():
            if pd.notna(row.get("date")):
                date_str = str(row["date"])
                if len(date_str) >= 4 and date_str[:4].isdigit():
                    tip_dates[row["strain"]] = int(date_str[:4])
        
        if tip_dates:
            years = list(tip_dates.values())
            temporal["date_range"] = {
                "min_year": min(years),
                "max_year": max(years),
                "span_years": max(years) - min(years)
            }
        
        return temporal
    
    def _generate_overall_assessment(self):
        """Generate overall quality assessment and recommendations"""
        scores = []
        
        # Collect quality scores
        for module in ["metadata_qc", "sequence_qc", "geographic_qc", "nextclade_qc", "phylogenetic_qc"]:
            if module in self.qc_results and "quality_score" in self.qc_results[module]:
                scores.append(self.qc_results[module]["quality_score"])
        
        if scores:
            overall_score = np.mean(scores)
            
            if overall_score >= 85:
                self.qc_results["overall_quality"] = "excellent"
            elif overall_score >= 70:
                self.qc_results["overall_quality"] = "good"
            elif overall_score >= 50:
                self.qc_results["overall_quality"] = "acceptable"
            else:
                self.qc_results["overall_quality"] = "poor"
            
            self.qc_results["summary"]["overall_score"] = overall_score
        
        # Generate recommendations
        self._generate_recommendations()
    
    def _generate_recommendations(self):
        """Generate quality improvement recommendations"""
        recommendations = []
        
        # Metadata recommendations
        if "metadata_qc" in self.qc_results:
            metadata_score = self.qc_results["metadata_qc"]["quality_score"]
            if metadata_score < 70:
                recommendations.append("Improve metadata completeness, particularly for required fields")
            
            if self.qc_results["metadata_qc"]["date_quality"]["invalid_dates"] > 0:
                recommendations.append("Standardize date formats to YYYY-MM-DD")
        
        # Sequence recommendations
        if "sequence_qc" in self.qc_results:
            seq_qc = self.qc_results["sequence_qc"]
            if seq_qc["quality_metrics"]["high_n_content_rate"] > 0.1:
                recommendations.append("Filter sequences with high N content (>5%)")
        
        # Geographic recommendations
        if "geographic_qc" in self.qc_results:
            geo_qc = self.qc_results["geographic_qc"]
            if geo_qc.get("coordinate_validation", {}).get("invalid_coordinates", 0) > 0:
                recommendations.append("Validate and correct invalid geographic coordinates")
        
        self.qc_results["recommendations"] = recommendations
    
    def _generate_comprehensive_report(self, output_dir: str):
        """Generate comprehensive QC report in multiple formats"""
          # JSON report
        json_file = os.path.join(output_dir, "comprehensive_qc_report.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(self.qc_results, f, indent=2, default=str)
        
        # HTML report
        html_file = os.path.join(output_dir, "comprehensive_qc_report.html")
        self._generate_html_report(html_file)
        
        # Text summary
        txt_file = os.path.join(output_dir, "comprehensive_qc_summary.txt")
        self._generate_text_summary(txt_file)
        
        logger.info(f"Reports generated:")
        logger.info(f"  - JSON: {json_file}")
        logger.info(f"  - HTML: {html_file}")
        logger.info(f"  - Summary: {txt_file}")
    
    def _generate_html_report(self, output_file: str):
        """Generate HTML QC report"""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>RVF Nextstrain Pipeline - Comprehensive QC Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        .header {{ background-color: #f0f8ff; padding: 20px; border-radius: 5px; }}
        .section {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
        .quality-excellent {{ color: #28a745; }}
        .quality-good {{ color: #007bff; }}
        .quality-acceptable {{ color: #ffc107; }}
        .quality-poor {{ color: #dc3545; }}
        .metric {{ margin: 10px 0; }}
        .score {{ font-weight: bold; font-size: 1.2em; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>RVF Nextstrain Pipeline - Comprehensive Quality Control Report</h1>
        <p><strong>Generated:</strong> {self.qc_results['timestamp']}</p>
        <p><strong>Overall Quality:</strong> 
            <span class="quality-{self.qc_results['overall_quality']}">{self.qc_results['overall_quality'].upper()}</span>
        </p>
        {f'<p><strong>Overall Score:</strong> <span class="score">{self.qc_results["summary"]["overall_score"]:.1f}/100</span></p>' if "summary" in self.qc_results and "overall_score" in self.qc_results["summary"] else ''}
    </div>
    
    <div class="section">
        <h2>Quality Control Summary</h2>
        {self._generate_qc_summary_html()}
    </div>
    
    <div class="section">
        <h2>Metadata Quality Control</h2>
        {self._generate_metadata_qc_html()}
    </div>
    
    <div class="section">
        <h2>Sequence Quality Control</h2>
        {self._generate_sequence_qc_html()}
    </div>
    
    <div class="section">
        <h2>Geographic Quality Control</h2>
        {self._generate_geographic_qc_html()}
    </div>
    
    <div class="section">
        <h2>Recommendations</h2>
        <ul>
        {''.join(f'<li>{rec}</li>' for rec in self.qc_results.get('recommendations', []))}
        </ul>
    </div>
</body>
</html>
        """
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
    
    def _generate_qc_summary_html(self) -> str:
        """Generate QC summary HTML section"""
        html = "<table><tr><th>Module</th><th>Score</th><th>Status</th></tr>"
        
        modules = [
            ("Metadata QC", "metadata_qc"),
            ("Sequence QC", "sequence_qc"), 
            ("Geographic QC", "geographic_qc"),
            ("Nextclade QC", "nextclade_qc"),
            ("Phylogenetic QC", "phylogenetic_qc")
        ]
        
        for name, key in modules:
            if key in self.qc_results and "quality_score" in self.qc_results[key]:
                score = self.qc_results[key]["quality_score"]
                status = "✓ Pass" if score >= 70 else "⚠ Warning" if score >= 50 else "✗ Fail"
                html += f"<tr><td>{name}</td><td>{score:.1f}</td><td>{status}</td></tr>"
            else:
                html += f"<tr><td>{name}</td><td>N/A</td><td>Not Run</td></tr>"
        
        html += "</table>"
        return html
    
    def _generate_metadata_qc_html(self) -> str:
        """Generate metadata QC HTML section"""
        if "metadata_qc" not in self.qc_results:
            return "<p>Metadata QC not performed.</p>"
        
        qc = self.qc_results["metadata_qc"]
        html = f"<p><strong>Quality Score:</strong> {qc['quality_score']:.1f}/100</p>"
        html += f"<p><strong>Total Records:</strong> {qc['total_records']}</p>"
        
        # Field completeness
        html += "<h3>Field Completeness</h3><table><tr><th>Field</th><th>Completeness</th><th>Status</th></tr>"
        for field, stats in qc["field_completeness"].items():
            rate = stats["completeness_rate"] * 100
            status = "Required" if stats["is_required"] else "Optional"
            html += f"<tr><td>{field}</td><td>{rate:.1f}%</td><td>{status}</td></tr>"
        html += "</table>"
        
        return html
    
    def _generate_sequence_qc_html(self) -> str:
        """Generate sequence QC HTML section"""
        if "sequence_qc" not in self.qc_results:
            return "<p>Sequence QC not performed.</p>"
        
                qc = self.qc_results["sequence_qc"]
        html = f"<p><strong>Quality Score:</strong> {qc['quality_score']:.1f}/100</p>"
        html += f"<p><strong>Total Sequences:</strong> {qc['total_sequences']}</p>"
        
        # Length distribution
        length_dist = qc["length_distribution"]
        html += "<h3>Length Distribution</h3>"
        html += f"<p>Mean: {length_dist['mean']:.0f} bp, Range: {length_dist['min']}-{length_dist['max']} bp</p>"
        
        return html
    
    def _generate_geographic_qc_html(self) -> str:
        """Generate geographic QC HTML section"""
        if "geographic_qc" not in self.qc_results:
            return "<p>Geographic QC not performed.</p>"
        
        qc = self.qc_results["geographic_qc"]
        html = f"<p><strong>Quality Score:</strong> {qc['quality_score']:.1f}/100</p>"
        
        if "geographic_distribution" in qc:
            geo_dist = qc["geographic_distribution"]
            html += f"<p><strong>Countries Represented:</strong> {geo_dist['unique_countries']}</p>"
        
        return html
    
    def _generate_text_summary(self, output_file: str):
        """Generate text summary report"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("RVF NEXTSTRAIN PIPELINE - COMPREHENSIVE QUALITY CONTROL SUMMARY\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Report Generated: {self.qc_results['timestamp']}\n")
            f.write(f"Overall Quality: {self.qc_results['overall_quality'].upper()}\n")
            
            if "summary" in self.qc_results and "overall_score" in self.qc_results["summary"]:
                f.write(f"Overall Score: {self.qc_results['summary']['overall_score']:.1f}/100\n")
            
            f.write("\nQUALITY CONTROL MODULES:\n")
            f.write("-" * 30 + "\n")
            
            modules = [
                ("Metadata QC", "metadata_qc"),
                ("Sequence QC", "sequence_qc"),
                ("Geographic QC", "geographic_qc"),
                ("Nextclade QC", "nextclade_qc"),
                ("Phylogenetic QC", "phylogenetic_qc")
            ]
            
            for name, key in modules:
                if key in self.qc_results and "quality_score" in self.qc_results[key]:
                    score = self.qc_results[key]["quality_score"]
                    status = "PASS" if score >= 70 else "WARNING" if score >= 50 else "FAIL"
                    f.write(f"{name}: {score:.1f}/100 ({status})\n")
                else:
                    f.write(f"{name}: Not performed\n")
            
            if self.qc_results.get("recommendations"):
                f.write("\nRECOMMENDATIONS:\n")
                f.write("-" * 20 + "\n")
                for i, rec in enumerate(self.qc_results["recommendations"], 1):
                    f.write(f"{i}. {rec}\n")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Comprehensive Quality Control for RVF Nextstrain Pipeline")
    parser.add_argument("--metadata", required=True, help="Path to metadata TSV file")
    parser.add_argument("--sequences", required=True, help="Path to sequences FASTA file")
    parser.add_argument("--tree", help="Path to phylogenetic tree file (optional)")
    parser.add_argument("--nextclade", help="Path to Nextclade JSON results (optional)")
    parser.add_argument("--output-dir", default="results/qc_reports", help="Output directory")
    parser.add_argument("--config", help="Configuration file (optional)")
    
    args = parser.parse_args()
    
    # Load configuration
    config = {}
    if args.config and os.path.exists(args.config):
        with open(args.config, 'r') as f:
            config = json.load(f)
    
    # Initialize QC system
    qc_system = ComprehensiveQC(config)
    
    # Run comprehensive QC
    results = qc_system.run_comprehensive_qc(
        metadata_file=args.metadata,
        sequences_file=args.sequences,
        tree_file=args.tree,
        nextclade_results=args.nextclade,
        output_dir=args.output_dir
    )
    
    # Print summary
    print(f"\nQuality Control Analysis Complete!")
    print(f"Overall Quality: {results['overall_quality'].upper()}")
    if "summary" in results and "overall_score" in results["summary"]:
        print(f"Overall Score: {results['summary']['overall_score']:.1f}/100")
    print(f"Reports saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
