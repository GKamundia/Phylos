#!/usr/bin/env python3
"""
RVF Nextstrain Pipeline Performance Assessment
Comprehensive analysis of pipeline efficiency and data processing metrics
"""

import pandas as pd
import json
import os
from datetime import datetime
from pathlib import Path

class PipelinePerformanceAssessment:
    """Assess RVF Nextstrain Pipeline Performance"""
    
    def __init__(self):
        self.assessment_results = {}
        self.base_dir = Path(".")
    
    def analyze_data_flow(self):
        """Analyze data flow through pipeline stages"""
        print("Analyzing data flow through pipeline stages...")
        
        # Raw data metrics
        raw_metadata_file = self.base_dir / "data/metadata/raw/rvf_metadata.tsv"
        raw_sequences_file = self.base_dir / "data/sequences/raw/rvf_sequences.fasta"
        
        raw_metadata_count = 0
        if raw_metadata_file.exists():
            raw_df = pd.read_csv(raw_metadata_file, sep='\t')
            raw_metadata_count = len(raw_df)
        
        raw_sequences_count = self._count_fasta_sequences(raw_sequences_file)
        
        # Enhanced metadata metrics
        enhanced_metadata_file = self.base_dir / "data/metadata/rvf_metadata.tsv"
        enhanced_metadata_count = 0
        enhanced_country_coverage = 0
        enhanced_coord_coverage = 0
        
        if enhanced_metadata_file.exists():
            enhanced_df = pd.read_csv(enhanced_metadata_file, sep='\t')
            enhanced_metadata_count = len(enhanced_df)
            if 'country' in enhanced_df.columns:
                enhanced_country_coverage = enhanced_df['country'].notna().sum()
            if 'latitude' in enhanced_df.columns and 'longitude' in enhanced_df.columns:
                enhanced_coord_coverage = (enhanced_df['latitude'].notna() & 
                                         enhanced_df['longitude'].notna()).sum()
        
        # Filtered data metrics
        filtered_metadata_file = self.base_dir / "results/filtered/rvf_metadata.tsv"
        filtered_sequences_file = self.base_dir / "results/filtered/rvf_filtered.fasta"
        
        filtered_metadata_count = 0
        filtered_sequences_count = 0
        
        if filtered_metadata_file.exists():
            filtered_df = pd.read_csv(filtered_metadata_file, sep='\t')
            filtered_metadata_count = len(filtered_df)
        
        filtered_sequences_count = self._count_fasta_sequences(filtered_sequences_file)
        
        # Calculate retention rates
        metadata_retention = (filtered_metadata_count / raw_metadata_count * 100) if raw_metadata_count > 0 else 0
        sequence_retention = (filtered_sequences_count / raw_sequences_count * 100) if raw_sequences_count > 0 else 0
        
        self.assessment_results['data_flow'] = {
            'raw_data': {
                'metadata_records': raw_metadata_count,
                'sequences': raw_sequences_count
            },
            'enhanced_data': {
                'metadata_records': enhanced_metadata_count,
                'country_coverage': enhanced_country_coverage,
                'coordinate_coverage': enhanced_coord_coverage,
                'country_extraction_rate': round(enhanced_country_coverage / enhanced_metadata_count * 100, 2) if enhanced_metadata_count > 0 else 0,
                'coordinate_mapping_rate': round(enhanced_coord_coverage / enhanced_metadata_count * 100, 2) if enhanced_metadata_count > 0 else 0
            },
            'filtered_data': {
                'metadata_records': filtered_metadata_count,
                'sequences': filtered_sequences_count,
                'metadata_retention_rate': round(metadata_retention, 2),
                'sequence_retention_rate': round(sequence_retention, 2)
            }
        }
    
    def analyze_processing_efficiency(self):
        """Analyze processing efficiency and bottlenecks"""
        print("Analyzing processing efficiency...")
        
        # Check log files for performance metrics
        log_dir = self.base_dir / "logs"
        benchmarks_dir = self.base_dir / "benchmarks"
        
        processing_stages = {}
        
        # Analyze benchmark files if available
        if benchmarks_dir.exists():
            for benchmark_file in benchmarks_dir.glob("*.txt"):
                stage_name = benchmark_file.stem
                try:
                    with open(benchmark_file, 'r') as f:
                        content = f.read()
                        # Extract timing information (simplified)
                        processing_stages[stage_name] = {
                            'benchmark_file': str(benchmark_file),
                            'has_performance_data': True
                        }
                except Exception as e:
                    processing_stages[stage_name] = {
                        'benchmark_file': str(benchmark_file),
                        'has_performance_data': False,
                        'error': str(e)
                    }
        
        # Check critical files and their sizes
        critical_files = {
            'raw_sequences': 'data/sequences/raw/rvf_sequences.fasta',
            'enhanced_metadata': 'data/metadata/rvf_metadata.tsv',
            'filtered_metadata': 'results/filtered/rvf_metadata.tsv',
            'filtered_sequences': 'results/filtered/rvf_filtered.fasta'
        }
        
        file_analysis = {}
        for name, filepath in critical_files.items():
            full_path = self.base_dir / filepath
            if full_path.exists():
                file_size = full_path.stat().st_size
                file_analysis[name] = {
                    'exists': True,
                    'size_bytes': file_size,
                    'size_mb': round(file_size / 1024 / 1024, 2),
                    'last_modified': datetime.fromtimestamp(full_path.stat().st_mtime).isoformat()
                }
            else:
                file_analysis[name] = {
                    'exists': False,
                    'size_bytes': 0,
                    'size_mb': 0
                }
        
        self.assessment_results['processing_efficiency'] = {
            'processing_stages': processing_stages,
            'file_analysis': file_analysis,
            'pipeline_completeness': self._assess_pipeline_completeness(file_analysis)
        }
    
    def analyze_quality_metrics(self):
        """Analyze quality control metrics"""
        print("Analyzing quality metrics...")
        
        # Load QC results if available
        qc_reports_dir = self.base_dir / "results/qc_reports"
        quality_assessment = {}
        
        if qc_reports_dir.exists():
            # Check for simple QC results
            simple_qc_file = qc_reports_dir / "simple_qc_report.json"
            if simple_qc_file.exists():
                try:
                    with open(simple_qc_file, 'r', encoding='utf-8') as f:
                        qc_data = json.load(f)
                        quality_assessment['simple_qc'] = {
                            'overall_quality': qc_data.get('overall_quality', 'unknown'),
                            'metadata_score': qc_data.get('metadata_qc', {}).get('quality_score', 0),
                            'sequence_score': qc_data.get('sequence_qc', {}).get('quality_score', 0),
                            'total_records': qc_data.get('metadata_qc', {}).get('total_records', 0),
                            'total_sequences': qc_data.get('sequence_qc', {}).get('total_sequences', 0)
                        }
                except Exception as e:
                    quality_assessment['simple_qc'] = {'error': str(e)}
        
        # Assess data completeness from the latest QC
        if 'simple_qc' in quality_assessment:
            qc_simple = quality_assessment['simple_qc']
            quality_assessment['overall_assessment'] = {
                'data_acquisition': 'excellent' if qc_simple.get('total_records', 0) > 1000 else 'good',
                'metadata_quality': self._score_to_quality(qc_simple.get('metadata_score', 0)),
                'sequence_quality': self._score_to_quality(qc_simple.get('sequence_score', 0)),
                'pipeline_status': 'partial' if qc_simple.get('total_sequences', 0) == 0 else 'complete'
            }
        
        self.assessment_results['quality_metrics'] = quality_assessment
    
    def generate_performance_report(self, output_dir):
        """Generate comprehensive performance report"""
        print("Generating performance report...")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Add metadata
        self.assessment_results['assessment_metadata'] = {
            'timestamp': datetime.now().isoformat(),
            'pipeline_version': 'Enhanced RVF Nextstrain Pipeline v2.0',
            'assessment_type': 'comprehensive_performance'
        }
          # JSON report
        json_file = os.path.join(output_dir, "pipeline_performance_assessment.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(self.assessment_results, f, indent=2, default=self._json_serializer)
        
        # Text summary
        summary_file = os.path.join(output_dir, "pipeline_performance_summary.txt")
        self._generate_text_summary(summary_file)
        
        print(f"Performance assessment reports generated in: {output_dir}")
    
    def _count_fasta_sequences(self, filepath):
        """Count sequences in FASTA file"""
        if not Path(filepath).exists():
            return 0
        
        try:
            with open(filepath, 'r') as f:
                count = sum(1 for line in f if line.startswith('>'))
            return count
        except Exception:
            return 0
    
    def _assess_pipeline_completeness(self, file_analysis):
        """Assess overall pipeline completeness"""
        required_files = ['raw_sequences', 'enhanced_metadata']
        optional_files = ['filtered_metadata', 'filtered_sequences']
        
        required_complete = all(file_analysis.get(f, {}).get('exists', False) for f in required_files)
        optional_complete = all(file_analysis.get(f, {}).get('exists', False) for f in optional_files)
        
        if required_complete and optional_complete:
            return 'complete'
        elif required_complete:
            return 'partial'
        else:
            return 'incomplete'
    
    def _score_to_quality(self, score):
        """Convert numeric score to quality label"""
        if score >= 85:
            return 'excellent'
        elif score >= 70:
            return 'good'
        elif score >= 50:
            return 'acceptable'
        else:
            return 'poor'
    
    def _json_serializer(self, obj):
        """JSON serializer for numpy types"""
        if hasattr(obj, 'item'):
            return obj.item()
        elif hasattr(obj, '__int__'):
            return int(obj)
        elif hasattr(obj, '__float__'):
            return float(obj)
        return str(obj)
    
    def _generate_text_summary(self, output_file):
        """Generate text summary of performance assessment"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("RVF NEXTSTRAIN PIPELINE - PERFORMANCE ASSESSMENT\n")
            f.write("=" * 55 + "\n\n")
            
            f.write(f"Assessment Date: {self.assessment_results['assessment_metadata']['timestamp']}\n")
            f.write(f"Pipeline Version: {self.assessment_results['assessment_metadata']['pipeline_version']}\n\n")
            
            # Data Flow Summary
            if 'data_flow' in self.assessment_results:
                data_flow = self.assessment_results['data_flow']
                f.write("DATA FLOW ANALYSIS\n")
                f.write("-" * 20 + "\n")
                
                raw = data_flow['raw_data']
                enhanced = data_flow['enhanced_data']
                filtered = data_flow['filtered_data']
                
                f.write(f"Raw Data Input:\n")
                f.write(f"  Metadata Records: {raw['metadata_records']:,}\n")
                f.write(f"  Sequences: {raw['sequences']:,}\n\n")
                
                f.write(f"Enhanced Processing:\n")
                f.write(f"  Country Extraction Rate: {enhanced['country_extraction_rate']}%\n")
                f.write(f"  Coordinate Mapping Rate: {enhanced['coordinate_mapping_rate']}%\n")
                f.write(f"  Records with Countries: {enhanced['country_coverage']:,}\n")
                f.write(f"  Records with Coordinates: {enhanced['coordinate_coverage']:,}\n\n")
                
                f.write(f"Filtered Output:\n")
                f.write(f"  Metadata Retention: {filtered['metadata_retention_rate']}%\n")
                f.write(f"  Sequence Retention: {filtered['sequence_retention_rate']}%\n")
                f.write(f"  Final Metadata Records: {filtered['metadata_records']:,}\n")
                f.write(f"  Final Sequences: {filtered['sequences']:,}\n\n")
            
            # Processing Efficiency Summary
            if 'processing_efficiency' in self.assessment_results:
                proc_eff = self.assessment_results['processing_efficiency']
                f.write("PROCESSING EFFICIENCY\n")
                f.write("-" * 20 + "\n")
                f.write(f"Pipeline Completeness: {proc_eff['pipeline_completeness'].upper()}\n")
                
                file_analysis = proc_eff['file_analysis']
                f.write("Critical Files Status:\n")
                for file_name, file_info in file_analysis.items():
                    status = "✓" if file_info['exists'] else "✗"
                    size = f"({file_info['size_mb']} MB)" if file_info['exists'] else ""
                    f.write(f"  {status} {file_name}: {size}\n")
                f.write("\n")
            
            # Quality Metrics Summary
            if 'quality_metrics' in self.assessment_results:
                quality = self.assessment_results['quality_metrics']
                f.write("QUALITY ASSESSMENT\n")
                f.write("-" * 18 + "\n")
                
                if 'overall_assessment' in quality:
                    overall = quality['overall_assessment']
                    f.write(f"Data Acquisition: {overall['data_acquisition'].upper()}\n")
                    f.write(f"Metadata Quality: {overall['metadata_quality'].upper()}\n")
                    f.write(f"Sequence Quality: {overall['sequence_quality'].upper()}\n")
                    f.write(f"Pipeline Status: {overall['pipeline_status'].upper()}\n")
                
                if 'simple_qc' in quality and 'error' not in quality['simple_qc']:
                    qc = quality['simple_qc']
                    f.write(f"Overall QC Rating: {qc['overall_quality'].upper()}\n")
                    f.write(f"Metadata Score: {qc['metadata_score']}/100\n")
                    f.write(f"Sequence Score: {qc['sequence_score']}/100\n")

def main():
    """Run pipeline performance assessment"""
    print("Starting pipeline performance assessment...")
    assessment = PipelinePerformanceAssessment()
    
    # Run all analyses
    assessment.analyze_data_flow()
    assessment.analyze_processing_efficiency()
    assessment.analyze_quality_metrics()
    
    # Generate reports
    assessment.generate_performance_report("results/performance_assessment")
    print("Pipeline performance assessment completed!")

if __name__ == "__main__":
    main()
