#!/usr/bin/env python3
"""
Simple Quality Control for RVF Nextstrain Pipeline
Focused on essential QC functionality without complex formatting
"""

import argparse
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
from datetime import datetime
from pathlib import Path

class SimpleQC:
    """Simple Quality Control for RVF Nextstrain Pipeline"""
    
    def __init__(self):
        self.qc_results = {}
    
    def analyze_metadata(self, metadata_file):
        """Analyze metadata quality"""
        print("Analyzing metadata quality...")
        
        df = pd.read_csv(metadata_file, sep='\t')
        total_records = len(df)
        
        # Field completeness analysis
        completeness = {}
        required_fields = ['strain', 'virus', 'accession', 'date', 'country']
        
        for field in required_fields:
            if field in df.columns:
                non_null = df[field].notna().sum()
                completeness[field] = {
                    'complete_count': int(non_null),
                    'completeness_rate': round(float(non_null / total_records * 100), 2)
                }
            else:
                completeness[field] = {
                    'complete_count': 0,
                    'completeness_rate': 0.0
                }
        
        # Date validation
        date_issues = 0
        if 'date' in df.columns:
            for date_val in df['date'].dropna():
                if not self._is_valid_date(str(date_val)):
                    date_issues += 1
        
        # Geographic distribution
        countries = df['country'].value_counts().to_dict() if 'country' in df.columns else {}
        
        self.qc_results['metadata_qc'] = {
            'total_records': total_records,
            'field_completeness': completeness,
            'date_issues': date_issues,
            'geographic_distribution': {
                'unique_countries': len(countries),
                'country_counts': dict(list(countries.items())[:10])  # Top 10
            },
            'quality_score': self._calculate_metadata_score(completeness, date_issues, total_records)
        }
    
    def analyze_sequences(self, sequences_file):
        """Analyze sequence quality"""
        print("Analyzing sequence quality...")
        
        sequences = list(SeqIO.parse(sequences_file, "fasta"))
        total_sequences = len(sequences)
        
        if total_sequences == 0:
            self.qc_results['sequence_qc'] = {
                'total_sequences': 0,
                'quality_score': 0
            }
            return
        
        # Length analysis
        lengths = [len(seq.seq) for seq in sequences]
        
        # GC content analysis
        gc_contents = []
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            if seq_str and not seq_str.replace('N', '').replace('-', '') == '':
                gc_contents.append(gc_fraction(seq_str) * 100)
        
        # N content analysis
        n_contents = []
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            if len(seq_str) > 0:
                n_content = seq_str.count('N') / len(seq_str) * 100
                n_contents.append(n_content)
        
        self.qc_results['sequence_qc'] = {
            'total_sequences': total_sequences,
            'length_distribution': {
                'mean': float(np.mean(lengths)) if lengths else 0,
                'median': float(np.median(lengths)) if lengths else 0,
                'min': int(min(lengths)) if lengths else 0,
                'max': int(max(lengths)) if lengths else 0,
                'std': float(np.std(lengths)) if lengths else 0
            },
            'gc_content': {
                'mean': float(np.mean(gc_contents)) if gc_contents else 0,
                'std': float(np.std(gc_contents)) if gc_contents else 0
            },
            'n_content': {
                'mean': float(np.mean(n_contents)) if n_contents else 0,
                'max': float(max(n_contents)) if n_contents else 0
            },
            'quality_score': self._calculate_sequence_score(lengths, gc_contents, n_contents)
        }
    
    def generate_reports(self, output_dir):
        """Generate QC reports"""
        print("Generating QC reports...")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Add timestamp and overall assessment
        self.qc_results['timestamp'] = datetime.now().isoformat()
        self.qc_results['overall_quality'] = self._assess_overall_quality()
        
        # JSON report
        json_file = os.path.join(output_dir, "simple_qc_report.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(self.qc_results, f, indent=2)
        
        # Text summary
        summary_file = os.path.join(output_dir, "simple_qc_summary.txt")
        self._generate_text_summary(summary_file)
        
        print(f"QC reports generated in: {output_dir}")
        print(f"Overall Quality Assessment: {self.qc_results['overall_quality']}")
    
    def _is_valid_date(self, date_str):
        """Validate date format"""
        valid_formats = ['%Y-%m-%d', '%Y-%m', '%Y']
        for fmt in valid_formats:
            try:
                datetime.strptime(date_str, fmt)
                return True
            except ValueError:
                continue
        return False
    
    def _calculate_metadata_score(self, completeness, date_issues, total_records):
        """Calculate metadata quality score"""
        if total_records == 0:
            return 0
        
        # Base score from field completeness
        total_completeness = sum(field['completeness_rate'] for field in completeness.values())
        avg_completeness = total_completeness / len(completeness)
        
        # Penalty for date issues
        date_penalty = (date_issues / total_records) * 100 if total_records > 0 else 0
        
        score = max(0, avg_completeness - date_penalty)
        return round(score, 1)
    
    def _calculate_sequence_score(self, lengths, gc_contents, n_contents):
        """Calculate sequence quality score"""
        if not lengths:
            return 0
        
        score = 100
        
        # Length consistency check
        if len(lengths) > 1:
            length_cv = np.std(lengths) / np.mean(lengths) if np.mean(lengths) > 0 else 1
            if length_cv > 0.5:  # High variation in lengths
                score -= 20
        
        # GC content check
        if gc_contents:
            avg_gc = np.mean(gc_contents)
            if avg_gc < 30 or avg_gc > 70:  # Unusual GC content
                score -= 15
        
        # N content check
        if n_contents:
            avg_n = np.mean(n_contents)
            if avg_n > 5:  # High N content
                score -= 25
        
        return round(max(0, score), 1)
    
    def _assess_overall_quality(self):
        """Assess overall quality"""
        scores = []
        
        if 'metadata_qc' in self.qc_results:
            scores.append(self.qc_results['metadata_qc']['quality_score'])
        
        if 'sequence_qc' in self.qc_results:
            scores.append(self.qc_results['sequence_qc']['quality_score'])
        
        if not scores:
            return "unknown"
        
        avg_score = np.mean(scores)
        
        if avg_score >= 85:
            return "excellent"
        elif avg_score >= 70:
            return "good"
        elif avg_score >= 50:
            return "acceptable"
        else:
            return "poor"
    
    def _generate_text_summary(self, output_file):
        """Generate text summary report"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("RVF NEXTSTRAIN PIPELINE - QUALITY CONTROL SUMMARY\n")
            f.write("=" * 55 + "\n\n")
            
            f.write(f"Report Generated: {self.qc_results['timestamp']}\n")
            f.write(f"Overall Quality: {self.qc_results['overall_quality'].upper()}\n\n")
            
            # Metadata QC Summary
            if 'metadata_qc' in self.qc_results:
                meta_qc = self.qc_results['metadata_qc']
                f.write("METADATA QUALITY CONTROL\n")
                f.write("-" * 25 + "\n")
                f.write(f"Total Records: {meta_qc['total_records']}\n")
                f.write(f"Quality Score: {meta_qc['quality_score']}/100\n")
                f.write(f"Date Issues: {meta_qc['date_issues']}\n")
                f.write(f"Countries Represented: {meta_qc['geographic_distribution']['unique_countries']}\n\n")
                
                f.write("Field Completeness:\n")
                for field, stats in meta_qc['field_completeness'].items():
                    f.write(f"  {field}: {stats['completeness_rate']}% ({stats['complete_count']} records)\n")
                f.write("\n")
            
            # Sequence QC Summary
            if 'sequence_qc' in self.qc_results:
                seq_qc = self.qc_results['sequence_qc']
                f.write("SEQUENCE QUALITY CONTROL\n")
                f.write("-" * 25 + "\n")
                f.write(f"Total Sequences: {seq_qc['total_sequences']}\n")
                f.write(f"Quality Score: {seq_qc['quality_score']}/100\n")
                
                if 'length_distribution' in seq_qc:
                    length_dist = seq_qc['length_distribution']
                    f.write(f"Length Range: {length_dist['min']}-{length_dist['max']} bp\n")
                    f.write(f"Mean Length: {length_dist['mean']:.0f} bp\n")
                
                if 'gc_content' in seq_qc:
                    gc = seq_qc['gc_content']
                    f.write(f"Average GC Content: {gc['mean']:.1f}%\n")
                
                if 'n_content' in seq_qc:
                    n_cont = seq_qc['n_content']
                    f.write(f"Average N Content: {n_cont['mean']:.1f}%\n")

def main():
    parser = argparse.ArgumentParser(description="Simple Quality Control for RVF Nextstrain Pipeline")
    parser.add_argument("--metadata", required=True, help="Path to metadata TSV file")
    parser.add_argument("--sequences", required=True, help="Path to sequences FASTA file")
    parser.add_argument("--output-dir", default="qc_reports", help="Output directory")
    
    args = parser.parse_args()
    
    # Create QC system
    qc_system = SimpleQC()
    
    # Run QC analysis
    qc_system.analyze_metadata(args.metadata)
    qc_system.analyze_sequences(args.sequences)
    
    # Generate reports
    qc_system.generate_reports(args.output_dir)

if __name__ == "__main__":
    main()
