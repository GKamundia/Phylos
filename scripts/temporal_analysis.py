#!/usr/bin/env python3
"""
Temporal analysis and molecular clock estimation for segmented viruses
Specialized for RVF and other segmented pathogens with enhanced temporal features
"""

import argparse
import logging
import os
import sys
import json
import subprocess
from pathlib import Path
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
from Bio.Phylo import TreeConstruction
import warnings
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class TemporalAnalysisEngine:
    """Comprehensive temporal analysis for phylogenetic data"""
    
    def __init__(self, output_dir=None):
        self.output_dir = output_dir or "temporal_analysis"
        self.results = {}
        self.plots_generated = []
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
    def analyze_temporal_patterns(self, metadata_file, tree_file=None):
        """
        Comprehensive temporal pattern analysis
        
        Args:
            metadata_file: Path to metadata with dates
            tree_file: Optional phylogenetic tree for advanced analysis
        """
        logger.info("Starting comprehensive temporal analysis")
        
        # Load and prepare data
        metadata = self._load_and_prepare_metadata(metadata_file)
        
        if metadata is None or len(metadata) == 0:
            logger.error("No valid metadata loaded")
            return None
            
        # Basic temporal statistics
        self._calculate_basic_temporal_stats(metadata)
        
        # Sampling pattern analysis
        self._analyze_sampling_patterns(metadata)
        
        # Geographic temporal patterns
        self._analyze_geographic_temporal_patterns(metadata)
        
        # Segment-specific temporal patterns (for segmented viruses)
        self._analyze_segment_temporal_patterns(metadata)
        
        # Generate visualizations
        self._generate_temporal_plots(metadata)
        
        # Advanced phylogenetic temporal analysis if tree provided
        if tree_file and os.path.exists(tree_file):
            self._advanced_phylogenetic_temporal_analysis(tree_file, metadata)
            
        # Save results
        self._save_results()
        
        return self.results
        
    def _load_and_prepare_metadata(self, metadata_file):
        """Load and prepare metadata for temporal analysis"""
        try:
            # Load metadata
            metadata = pd.read_csv(metadata_file, sep='\t')
            logger.info(f"Loaded {len(metadata)} records from metadata")
              # Prepare date column - check for various date column names
            date_columns = ['date', 'Date', 'Release_Date', 'collection_date', 'Collection_Date']
            date_col = None
            
            for col in date_columns:
                if col in metadata.columns:
                    date_col = col
                    break
                    
            if date_col is None:
                logger.error(f"No date column found in metadata. Available columns: {list(metadata.columns)}")
                return None
                
            logger.info(f"Using date column: {date_col}")
              # Rename to standard 'date' column
            if date_col != 'date':
                metadata['date'] = metadata[date_col]
                
            # Convert dates with multiple format support - try pandas auto parsing first
            try:
                metadata['date'] = pd.to_datetime(metadata['date'], errors='coerce')
                logger.info("Successfully used pandas automatic date parsing")
            except Exception:
                # Fallback to specific format parsing
                date_formats = ['%Y-%m-%d', '%Y-%m', '%Y', '%d/%m/%Y', '%m/%d/%Y', '%d-%b-%Y', '%d-%B-%Y']
                metadata['date'] = None
                
                for fmt in date_formats:
                    mask = metadata['date'].isna()
                    if mask.any():
                        try:
                            parsed_dates = pd.to_datetime(metadata.loc[mask, date_col], format=fmt, errors='coerce')
                            metadata.loc[mask, 'date'] = parsed_dates
                            if not metadata['date'].isna().all():
                                logger.info(f"Successfully parsed dates using format: {fmt}")
                        except Exception:
                            continue
                
            # Remove records without valid dates
            initial_count = len(metadata)
            metadata = metadata.dropna(subset=['date'])
            final_count = len(metadata)
            
            logger.info(f"Successfully parsed dates for {final_count}/{initial_count} records")
            
            # Add derived temporal columns
            metadata['year'] = metadata['date'].dt.year
            metadata['month'] = metadata['date'].dt.month
            metadata['day_of_year'] = metadata['date'].dt.dayofyear
            metadata['quarter'] = metadata['date'].dt.quarter
            
            return metadata
            
        except Exception as e:
            logger.error(f"Failed to load metadata: {e}")
            return None
            
    def _calculate_basic_temporal_stats(self, metadata):
        """Calculate basic temporal statistics"""
        logger.info("Calculating basic temporal statistics")
        
        date_range = metadata['date'].max() - metadata['date'].min()
        
        basic_stats = {
            'total_sequences': len(metadata),
            'date_range_days': date_range.days,
            'date_range_years': date_range.days / 365.25,
            'earliest_date': metadata['date'].min().isoformat(),
            'latest_date': metadata['date'].max().isoformat(),
            'year_range': f"{metadata['year'].min()}-{metadata['year'].max()}",
            'unique_years': metadata['year'].nunique(),
            'unique_months': len(metadata[['year', 'month']].drop_duplicates())
        }
        
        self.results['basic_temporal_stats'] = basic_stats
        logger.info(f"Temporal range: {basic_stats['date_range_years']:.1f} years")
        
    def _analyze_sampling_patterns(self, metadata):
        """Analyze temporal sampling patterns"""
        logger.info("Analyzing sampling patterns")
        
        # Yearly sampling
        yearly_counts = metadata['year'].value_counts().sort_index()
        
        # Monthly sampling
        monthly_counts = metadata.groupby(['year', 'month']).size()
        
        # Quarterly sampling
        quarterly_counts = metadata.groupby(['year', 'quarter']).size()
        
        sampling_patterns = {
            'yearly_sampling': yearly_counts.to_dict(),
            'peak_year': yearly_counts.idxmax(),
            'peak_year_count': yearly_counts.max(),
            'avg_sequences_per_year': yearly_counts.mean(),
            'median_sequences_per_year': yearly_counts.median(),
            'sampling_variance': yearly_counts.var(),
            'years_with_sampling': yearly_counts[yearly_counts > 0].index.tolist(),
            'gaps_in_sampling': self._identify_sampling_gaps(yearly_counts)
        }
        
        # Seasonal patterns
        seasonal_patterns = metadata.groupby(metadata['date'].dt.month).size()
        sampling_patterns['seasonal_distribution'] = seasonal_patterns.to_dict()
        sampling_patterns['peak_month'] = seasonal_patterns.idxmax()
        
        self.results['sampling_patterns'] = sampling_patterns
        
    def _identify_sampling_gaps(self, yearly_counts):
        """Identify gaps in temporal sampling"""
        years = yearly_counts.index
        year_range = range(years.min(), years.max() + 1)
        missing_years = [year for year in year_range if year not in years]
        
        # Identify consecutive gaps
        gaps = []
        if missing_years:
            current_gap = [missing_years[0]]
            
            for i in range(1, len(missing_years)):
                if missing_years[i] == missing_years[i-1] + 1:
                    current_gap.append(missing_years[i])
                else:
                    gaps.append(current_gap)
                    current_gap = [missing_years[i]]
            gaps.append(current_gap)
            
        return gaps
        
    def _analyze_geographic_temporal_patterns(self, metadata):
        """Analyze temporal patterns by geographic location"""
        logger.info("Analyzing geographic temporal patterns")
        
        geographic_temporal = {}
        
        # Analyze by country if available
        if 'country' in metadata.columns:
            country_temporal = {}
            
            for country in metadata['country'].unique():
                if pd.isna(country):
                    continue
                    
                country_data = metadata[metadata['country'] == country]
                country_temporal[country] = {
                    'count': len(country_data),
                    'date_range': (country_data['date'].max() - country_data['date'].min()).days,
                    'first_sample': country_data['date'].min().isoformat(),
                    'last_sample': country_data['date'].max().isoformat(),
                    'years_active': country_data['year'].nunique(),
                    'yearly_distribution': country_data['year'].value_counts().to_dict()
                }
                
            geographic_temporal['by_country'] = country_temporal
            
        # Analyze by division if available
        if 'division' in metadata.columns:
            division_temporal = {}
            
            for division in metadata['division'].unique():
                if pd.isna(division):
                    continue
                    
                division_data = metadata[metadata['division'] == division]
                if len(division_data) > 1:  # Only include divisions with multiple samples
                    division_temporal[division] = {
                        'count': len(division_data),
                        'date_range': (division_data['date'].max() - division_data['date'].min()).days,
                        'years_active': division_data['year'].nunique()
                    }
                    
            geographic_temporal['by_division'] = division_temporal
            
        self.results['geographic_temporal_patterns'] = geographic_temporal
        
    def _analyze_segment_temporal_patterns(self, metadata):
        """Analyze temporal patterns by genomic segment (for segmented viruses)"""
        logger.info("Analyzing segment-specific temporal patterns")
        
        if 'segment' not in metadata.columns:
            logger.info("No segment information available, skipping segment analysis")
            return
            
        segment_temporal = {}
        
        for segment in metadata['segment'].unique():
            if pd.isna(segment):
                continue
                
            segment_data = metadata[metadata['segment'] == segment]
            
            segment_temporal[segment] = {
                'count': len(segment_data),
                'date_range': (segment_data['date'].max() - segment_data['date'].min()).days,
                'first_sample': segment_data['date'].min().isoformat(),
                'last_sample': segment_data['date'].max().isoformat(),
                'years_active': segment_data['year'].nunique(),
                'yearly_distribution': segment_data['year'].value_counts().to_dict(),
                'countries': segment_data['country'].nunique() if 'country' in segment_data.columns else 0
            }
            
        self.results['segment_temporal_patterns'] = segment_temporal
        
    def _generate_temporal_plots(self, metadata):
        """Generate temporal visualization plots"""
        logger.info("Generating temporal plots")
        
        # Set plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 1. Yearly sampling distribution
        self._plot_yearly_distribution(metadata)
        
        # 2. Monthly sampling patterns
        self._plot_monthly_patterns(metadata)
        
        # 3. Geographic temporal patterns
        self._plot_geographic_temporal(metadata)
        
        # 4. Segment temporal patterns (if applicable)
        if 'segment' in metadata.columns:
            self._plot_segment_temporal(metadata)
            
    def _plot_yearly_distribution(self, metadata):
        """Plot yearly sampling distribution"""
        try:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
            
            # Yearly counts
            yearly_counts = metadata['year'].value_counts().sort_index()
            yearly_counts.plot(kind='bar', ax=ax1, color='steelblue', alpha=0.7)
            ax1.set_title('Sequences by Year', fontsize=14, fontweight='bold')
            ax1.set_xlabel('Year')
            ax1.set_ylabel('Number of Sequences')
            ax1.tick_params(axis='x', rotation=45)
            
            # Cumulative sequences over time
            cumulative = yearly_counts.cumsum()
            cumulative.plot(kind='line', ax=ax2, color='darkgreen', linewidth=2, marker='o')
            ax2.set_title('Cumulative Sequences Over Time', fontsize=14, fontweight='bold')
            ax2.set_xlabel('Year')
            ax2.set_ylabel('Cumulative Sequences')
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plot_file = os.path.join(self.output_dir, 'yearly_distribution.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.plots_generated.append(plot_file)
            
        except Exception as e:
            logger.warning(f"Failed to generate yearly distribution plot: {e}")
            
    def _plot_monthly_patterns(self, metadata):
        """Plot monthly sampling patterns"""
        try:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Monthly distribution across all years
            monthly_dist = metadata['month'].value_counts().sort_index()
            month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
            
            monthly_dist.plot(kind='bar', ax=ax1, color='coral', alpha=0.7)
            ax1.set_title('Sampling Distribution by Month', fontsize=14, fontweight='bold')
            ax1.set_xlabel('Month')
            ax1.set_ylabel('Number of Sequences')
            ax1.set_xticklabels(month_names, rotation=45)
            
            # Heatmap of year vs month
            year_month_pivot = metadata.groupby(['year', 'month']).size().unstack(fill_value=0)
            sns.heatmap(year_month_pivot, annot=True, fmt='d', cmap='YlOrRd', ax=ax2)
            ax2.set_title('Sampling Intensity Heatmap (Year vs Month)', fontsize=14, fontweight='bold')
            ax2.set_xlabel('Month')
            ax2.set_ylabel('Year')
            
            plt.tight_layout()
            plot_file = os.path.join(self.output_dir, 'monthly_patterns.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.plots_generated.append(plot_file)
            
        except Exception as e:
            logger.warning(f"Failed to generate monthly patterns plot: {e}")
            
    def _plot_geographic_temporal(self, metadata):
        """Plot geographic temporal patterns"""
        try:
            if 'country' not in metadata.columns:
                return
                
            # Top countries by sequence count
            top_countries = metadata['country'].value_counts().head(10)
            
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Create temporal plot for top countries
            for country in top_countries.index:
                if pd.isna(country):
                    continue
                country_data = metadata[metadata['country'] == country]
                yearly_counts = country_data['year'].value_counts().sort_index()
                ax.plot(yearly_counts.index, yearly_counts.values, 
                       marker='o', linewidth=2, label=country, alpha=0.8)
                
            ax.set_title('Temporal Sampling Patterns by Country (Top 10)', fontsize=14, fontweight='bold')
            ax.set_xlabel('Year')
            ax.set_ylabel('Number of Sequences')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plot_file = os.path.join(self.output_dir, 'geographic_temporal.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.plots_generated.append(plot_file)
            
        except Exception as e:
            logger.warning(f"Failed to generate geographic temporal plot: {e}")
            
    def _plot_segment_temporal(self, metadata):
        """Plot segment-specific temporal patterns"""
        try:
            segments = metadata['segment'].unique()
            segments = [s for s in segments if not pd.isna(s)]
            
            if len(segments) <= 1:
                return
                
            fig, axes = plt.subplots(len(segments), 1, figsize=(12, 4*len(segments)))
            if len(segments) == 1:
                axes = [axes]
                
            for i, segment in enumerate(segments):
                segment_data = metadata[metadata['segment'] == segment]
                yearly_counts = segment_data['year'].value_counts().sort_index()
                
                yearly_counts.plot(kind='bar', ax=axes[i], color=f'C{i}', alpha=0.7)
                axes[i].set_title(f'Segment {segment} - Temporal Distribution', 
                                fontsize=12, fontweight='bold')
                axes[i].set_xlabel('Year')
                axes[i].set_ylabel('Number of Sequences')
                axes[i].tick_params(axis='x', rotation=45)
                
            plt.tight_layout()
            plot_file = os.path.join(self.output_dir, 'segment_temporal.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.plots_generated.append(plot_file)
            
        except Exception as e:
            logger.warning(f"Failed to generate segment temporal plot: {e}")
            
    def _advanced_phylogenetic_temporal_analysis(self, tree_file, metadata):
        """Perform advanced phylogenetic temporal analysis"""
        logger.info("Performing advanced phylogenetic temporal analysis")
        
        try:
            # Load tree
            tree = Phylo.read(tree_file, "newick")
            
            # Create strain to date mapping
            strain_dates = {}
            for _, row in metadata.iterrows():
                strain_dates[row['strain']] = row['date']
                
            # Calculate temporal signal statistics
            temporal_signal = self._calculate_temporal_signal(tree, strain_dates)
            
            self.results['phylogenetic_temporal_analysis'] = temporal_signal
            
        except Exception as e:
            logger.warning(f"Advanced phylogenetic temporal analysis failed: {e}")
            
    def _calculate_temporal_signal(self, tree, strain_dates):
        """Calculate temporal signal in phylogenetic tree"""
        try:
            # Get terminal nodes with dates
            terminals_with_dates = []
            
            for terminal in tree.get_terminals():
                if terminal.name in strain_dates:
                    terminals_with_dates.append({
                        'name': terminal.name,
                        'date': strain_dates[terminal.name],
                        'distance_to_root': tree.distance(terminal)
                    })
                    
            if len(terminals_with_dates) < 3:
                return {'error': 'Insufficient terminals with dates for analysis'}
                
            # Convert to DataFrame for analysis
            df = pd.DataFrame(terminals_with_dates)
            df['numeric_date'] = (df['date'] - df['date'].min()).dt.days
            
            # Calculate correlation between genetic distance and time
            correlation = df['distance_to_root'].corr(df['numeric_date'])
            
            return {
                'temporal_signal_correlation': correlation,
                'terminals_analyzed': len(terminals_with_dates),
                'date_range_analyzed': (df['date'].max() - df['date'].min()).days
            }
            
        except Exception as e:
            return {'error': f'Temporal signal calculation failed: {e}'}
            
    def _save_results(self):
        """Save analysis results to JSON file"""
        # Add metadata about the analysis
        self.results['analysis_metadata'] = {
            'timestamp': datetime.now().isoformat(),
            'output_directory': self.output_dir,
            'plots_generated': self.plots_generated
        }
        
        # Save to JSON
        results_file = os.path.join(self.output_dir, 'temporal_analysis_results.json')
        with open(results_file, 'w') as f:
            json.dump(self.results, f, indent=2, default=str)
            
        logger.info(f"Results saved to {results_file}")
        
        # Create summary report
        self._create_summary_report()
        
    def _create_summary_report(self):
        """Create a human-readable summary report"""
        report_file = os.path.join(self.output_dir, 'temporal_analysis_summary.txt')
        
        with open(report_file, 'w') as f:
            f.write("TEMPORAL ANALYSIS SUMMARY REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            # Basic statistics
            if 'basic_temporal_stats' in self.results:
                stats = self.results['basic_temporal_stats']
                f.write("BASIC TEMPORAL STATISTICS\n")
                f.write("-" * 30 + "\n")
                f.write(f"Total sequences: {stats.get('total_sequences', 'N/A')}\n")
                f.write(f"Date range: {stats.get('year_range', 'N/A')}\n")
                f.write(f"Duration: {stats.get('date_range_years', 'N/A'):.1f} years\n")
                f.write(f"Unique years: {stats.get('unique_years', 'N/A')}\n\n")
                
            # Sampling patterns
            if 'sampling_patterns' in self.results:
                patterns = self.results['sampling_patterns']
                f.write("SAMPLING PATTERNS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Peak sampling year: {patterns.get('peak_year', 'N/A')}\n")
                f.write(f"Peak year count: {patterns.get('peak_year_count', 'N/A')}\n")
                f.write(f"Average per year: {patterns.get('avg_sequences_per_year', 'N/A'):.1f}\n")
                f.write(f"Peak month: {patterns.get('peak_month', 'N/A')}\n\n")
                
            # Geographic patterns
            if 'geographic_temporal_patterns' in self.results:
                geo = self.results['geographic_temporal_patterns']
                if 'by_country' in geo:
                    f.write("GEOGRAPHIC DISTRIBUTION\n")
                    f.write("-" * 25 + "\n")
                    countries = list(geo['by_country'].keys())[:5]
                    f.write(f"Top countries: {', '.join(countries)}\n\n")
                    
            # Segment patterns
            if 'segment_temporal_patterns' in self.results:
                segments = self.results['segment_temporal_patterns']
                f.write("SEGMENT ANALYSIS\n")
                f.write("-" * 18 + "\n")
                for segment, data in segments.items():
                    f.write(f"Segment {segment}: {data.get('count', 'N/A')} sequences\n")
                f.write("\n")
                
            f.write(f"Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Plots generated: {len(self.plots_generated)}\n")
            
        logger.info(f"Summary report saved to {report_file}")


def main():
    parser = argparse.ArgumentParser(description="Comprehensive temporal analysis for phylogenetic data")
    parser.add_argument("metadata", help="Metadata file with temporal information")
    parser.add_argument("--tree", help="Phylogenetic tree file for advanced analysis")
    parser.add_argument("--output-dir", default="temporal_analysis", 
                       help="Output directory for results and plots")
    parser.add_argument("--no-plots", action='store_true',
                       help="Skip plot generation")
    
    args = parser.parse_args()
    
    try:
        # Initialize temporal analysis engine
        analyzer = TemporalAnalysisEngine(output_dir=args.output_dir)
        
        # Run comprehensive analysis
        results = analyzer.analyze_temporal_patterns(args.metadata, args.tree)
        
        if results:
            logger.info(f"Temporal analysis completed successfully")
            logger.info(f"Results saved to: {args.output_dir}")
        else:
            logger.error("Temporal analysis failed")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Temporal analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
