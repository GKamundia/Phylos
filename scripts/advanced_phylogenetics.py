#!/usr/bin/env python3
"""
Advanced phylogenetic tree building and temporal analysis for RVF Nextstrain
Implements multiple tree building algorithms and temporal analysis features
"""

import argparse
import logging
import os
import sys
import tempfile
import json
import subprocess
from pathlib import Path
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import ParsimonyScorer, ParsimonyTreeConstructor
import dendropy
from tqdm import tqdm

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AdvancedPhylogeneticAnalyzer:
    """Advanced phylogenetic analysis with multiple algorithms and temporal features"""
    
    def __init__(self, method='iqtree', threads=4, temp_dir=None):
        self.method = method.lower()
        self.threads = threads
        self.temp_dir = temp_dir or tempfile.gettempdir()
        self.tree_stats = {}
        self.temporal_stats = {}
        
    def build_tree(self, alignment_file, output_tree, metadata_file=None, model='auto'):
        """
        Build phylogenetic tree using specified method
        
        Args:
            alignment_file: Path to aligned sequences
            output_tree: Path to output tree file
            metadata_file: Optional metadata for temporal analysis
            model: Evolutionary model (auto, GTR, HKY, etc.)
        """
        logger.info(f"Building phylogenetic tree with {self.method}")
        
        if self.method == 'iqtree':
            self._build_tree_iqtree(alignment_file, output_tree, model)
        elif self.method == 'fasttree':
            self._build_tree_fasttree(alignment_file, output_tree)
        elif self.method == 'raxml':
            self._build_tree_raxml(alignment_file, output_tree, model)
        elif self.method == 'neighbor_joining':
            self._build_tree_nj(alignment_file, output_tree)
        elif self.method == 'maximum_parsimony':
            self._build_tree_mp(alignment_file, output_tree)
        elif self.method == 'auto':
            self._auto_select_tree_method(alignment_file, output_tree, model)
        else:
            raise ValueError(f"Unsupported tree building method: {self.method}")
            
        # Calculate tree statistics
        self._calculate_tree_stats(output_tree)
        
        # Perform temporal analysis if metadata provided
        if metadata_file:
            self._temporal_analysis(output_tree, metadata_file)
            
    def _build_tree_iqtree(self, alignment_file, output_tree, model):
        """Build tree using IQ-TREE"""
        try:
            prefix = os.path.splitext(output_tree)[0]
            
            cmd = [
                'iqtree2',
                '-s', alignment_file,
                '-pre', prefix,
                '-nt', str(self.threads),
                '-redo'
            ]
            
            if model != 'auto':
                cmd.extend(['-m', model])
            else:
                cmd.extend(['-m', 'MFP'])  # Model finder plus
                
            logger.info(f"Running IQ-TREE: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"IQ-TREE failed: {result.stderr}")
                raise RuntimeError(f"IQ-TREE failed: {result.stderr}")
                
            # Copy the resulting tree
            iqtree_output = f"{prefix}.treefile"
            if os.path.exists(iqtree_output):
                os.rename(iqtree_output, output_tree)
            else:
                raise RuntimeError("IQ-TREE did not produce expected output file")
                
        except FileNotFoundError:
            logger.warning("IQ-TREE not found, falling back to FastTree")
            self._build_tree_fasttree(alignment_file, output_tree)
            
    def _build_tree_fasttree(self, alignment_file, output_tree):
        """Build tree using FastTree"""
        try:
            cmd = ['fasttree', '-nt', '-gtr', alignment_file]
            logger.info(f"Running FastTree: {' '.join(cmd)}")
            
            with open(output_tree, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
                
            if result.returncode != 0:
                logger.error(f"FastTree failed: {result.stderr}")
                raise RuntimeError(f"FastTree failed: {result.stderr}")
                
        except FileNotFoundError:
            logger.warning("FastTree not found, falling back to neighbor joining")
            self._build_tree_nj(alignment_file, output_tree)
            
    def _build_tree_raxml(self, alignment_file, output_tree, model):
        """Build tree using RAxML"""
        try:
            prefix = os.path.splitext(output_tree)[0]
            working_dir = os.path.dirname(output_tree)
            
            cmd = [
                'raxml-ng',
                '--msa', alignment_file,
                '--model', model if model != 'auto' else 'GTR+G',
                '--prefix', prefix,
                '--threads', str(self.threads),
                '--force'
            ]
            
            logger.info(f"Running RAxML-NG: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, cwd=working_dir, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"RAxML-NG failed: {result.stderr}")
                raise RuntimeError(f"RAxML-NG failed: {result.stderr}")
                
            # Copy the resulting tree
            raxml_output = f"{prefix}.raxml.bestTree"
            if os.path.exists(raxml_output):
                os.rename(raxml_output, output_tree)
            else:
                raise RuntimeError("RAxML-NG did not produce expected output file")
                
        except FileNotFoundError:
            logger.warning("RAxML-NG not found, falling back to FastTree")
            self._build_tree_fasttree(alignment_file, output_tree)
            
    def _build_tree_nj(self, alignment_file, output_tree):
        """Build tree using Neighbor Joining (BioPython implementation)"""
        logger.info("Building tree with Neighbor Joining")
        
        from Bio import AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        
        # Load alignment
        alignment = AlignIO.read(alignment_file, "fasta")
        
        # Calculate distance matrix
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
        
        # Build tree
        constructor = DistanceTreeConstructor(calculator, 'nj')
        tree = constructor.build_tree(alignment)
        
        # Save tree
        Phylo.write(tree, output_tree, "newick")
        
    def _build_tree_mp(self, alignment_file, output_tree):
        """Build tree using Maximum Parsimony (BioPython implementation)"""
        logger.info("Building tree with Maximum Parsimony")
        
        from Bio import AlignIO
        from Bio.Phylo.TreeConstruction import ParsimonyScorer, ParsimonyTreeConstructor
        
        # Load alignment
        alignment = AlignIO.read(alignment_file, "fasta")
        
        # Build tree using parsimony
        scorer = ParsimonyScorer()
        constructor = ParsimonyTreeConstructor(scorer, 'nni')
        tree = constructor.build_tree(alignment)
        
        # Save tree
        Phylo.write(tree, output_tree, "newick")
        
    def _auto_select_tree_method(self, alignment_file, output_tree, model):
        """Automatically select best tree building method based on data characteristics"""
        from Bio import AlignIO
        
        # Analyze alignment characteristics
        alignment = AlignIO.read(alignment_file, "fasta")
        num_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        logger.info(f"Auto-selecting tree method for {num_sequences} sequences, length {alignment_length}")
        
        # Decision logic for method selection
        if num_sequences < 50:
            selected_method = 'iqtree'
        elif num_sequences < 200:
            selected_method = 'fasttree'
        else:
            # For very large datasets, use fastest methods
            selected_method = 'neighbor_joining'
            
        logger.info(f"Selected tree building method: {selected_method}")
        
        # Update method and build tree
        original_method = self.method
        self.method = selected_method
        self.build_tree(alignment_file, output_tree, model=model)
        self.method = original_method
        
    def _calculate_tree_stats(self, tree_file):
        """Calculate and store tree statistics"""
        try:
            tree = Phylo.read(tree_file, "newick")
            
            # Basic tree statistics
            terminals = tree.get_terminals()
            nonterminals = tree.get_nonterminals()
            
            self.tree_stats = {
                'num_terminals': len(terminals),
                'num_internal_nodes': len(nonterminals),
                'total_nodes': len(terminals) + len(nonterminals),
                'tree_depth': tree.depth(),
                'timestamp': datetime.now().isoformat()
            }
            
            # Calculate branch length statistics
            branch_lengths = []
            for node in tree.find_clades():
                if node.branch_length:
                    branch_lengths.append(node.branch_length)
                    
            if branch_lengths:
                self.tree_stats.update({
                    'total_branch_length': sum(branch_lengths),
                    'avg_branch_length': np.mean(branch_lengths),
                    'median_branch_length': np.median(branch_lengths),
                    'max_branch_length': max(branch_lengths),
                    'min_branch_length': min(branch_lengths)
                })
                
            logger.info(f"Tree statistics: {self.tree_stats}")
            
        except Exception as e:
            logger.warning(f"Could not calculate tree statistics: {e}")
            
    def _temporal_analysis(self, tree_file, metadata_file):
        """Perform temporal analysis on the phylogenetic tree"""
        logger.info("Performing temporal analysis")
        
        try:
            # Load metadata
            metadata = pd.read_csv(metadata_file, sep='\t')
            
            # Convert dates
            if 'date' in metadata.columns:
                metadata['date'] = pd.to_datetime(metadata['date'], errors='coerce')
                
                # Calculate temporal statistics
                date_range = metadata['date'].max() - metadata['date'].min()
                
                self.temporal_stats = {
                    'date_range_days': date_range.days,
                    'earliest_date': metadata['date'].min().isoformat() if not pd.isna(metadata['date'].min()) else None,
                    'latest_date': metadata['date'].max().isoformat() if not pd.isna(metadata['date'].max()) else None,
                    'num_dated_sequences': metadata['date'].notna().sum(),
                    'temporal_coverage': metadata['date'].notna().sum() / len(metadata)
                }
                
                # Calculate sampling frequency by year
                if not metadata['date'].isna().all():
                    metadata['year'] = metadata['date'].dt.year
                    year_counts = metadata['year'].value_counts().sort_index()
                    
                    self.temporal_stats['sampling_by_year'] = year_counts.to_dict()
                    self.temporal_stats['peak_sampling_year'] = year_counts.idxmax()
                    self.temporal_stats['avg_sequences_per_year'] = year_counts.mean()
                    
            logger.info(f"Temporal analysis complete: {self.temporal_stats}")
            
        except Exception as e:
            logger.warning(f"Temporal analysis failed: {e}")


class TemporalPhylogeneticAnalyzer:
    """Specialized analyzer for time-resolved phylogenetic analysis"""
    
    def __init__(self, threads=4):
        self.threads = threads
        self.clock_analysis = {}
        
    def estimate_molecular_clock(self, tree_file, metadata_file, output_dir):
        """
        Estimate molecular clock and perform time-resolved analysis
        
        Args:
            tree_file: Path to phylogenetic tree
            metadata_file: Path to metadata with dates
            output_dir: Directory for output files
        """
        logger.info("Estimating molecular clock")
        
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # Use TreeTime for molecular clock estimation
            self._run_treetime(tree_file, metadata_file, output_dir)
            
        except Exception as e:
            logger.warning(f"TreeTime analysis failed: {e}")
            # Fallback to basic temporal analysis
            self._basic_temporal_analysis(tree_file, metadata_file, output_dir)
            
    def _run_treetime(self, tree_file, metadata_file, output_dir):
        """Run TreeTime for molecular clock estimation"""
        try:
            cmd = [
                'treetime',
                '--tree', tree_file,
                '--dates', metadata_file,
                '--outdir', output_dir,
                '--clock-filter', '3'
            ]
            
            logger.info(f"Running TreeTime: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"TreeTime failed: {result.stderr}")
                raise RuntimeError(f"TreeTime failed: {result.stderr}")
                
            # Parse TreeTime results
            self._parse_treetime_results(output_dir)
            
        except FileNotFoundError:
            logger.warning("TreeTime not found, using basic temporal analysis")
            raise
            
    def _parse_treetime_results(self, output_dir):
        """Parse TreeTime output and extract clock statistics"""
        try:
            # Look for TreeTime output files
            results_files = [
                'molecular_clock.txt',
                'dates.tsv',
                'timetree.nexus'
            ]
            
            for file in results_files:
                filepath = os.path.join(output_dir, file)
                if os.path.exists(filepath):
                    logger.info(f"Found TreeTime output: {file}")
                    
            # TODO: Parse specific TreeTime results
            self.clock_analysis['treetime_completed'] = True
            self.clock_analysis['output_dir'] = output_dir
            
        except Exception as e:
            logger.warning(f"Could not parse TreeTime results: {e}")
            
    def _basic_temporal_analysis(self, tree_file, metadata_file, output_dir):
        """Perform basic temporal analysis without TreeTime"""
        logger.info("Performing basic temporal analysis")
        
        # Load metadata
        metadata = pd.read_csv(metadata_file, sep='\t')
        
        # Basic date analysis
        if 'date' in metadata.columns:
            metadata['date'] = pd.to_datetime(metadata['date'], errors='coerce')
            
            # Calculate basic temporal statistics
            self.clock_analysis = {
                'method': 'basic',
                'date_range': (metadata['date'].max() - metadata['date'].min()).days,
                'num_sequences': len(metadata),
                'num_dated': metadata['date'].notna().sum(),
                'earliest_sample': metadata['date'].min().isoformat() if not pd.isna(metadata['date'].min()) else None,
                'latest_sample': metadata['date'].max().isoformat() if not pd.isna(metadata['date'].max()) else None
            }
            
            # Save basic results
            output_file = os.path.join(output_dir, 'temporal_analysis.json')
            with open(output_file, 'w') as f:
                json.dump(self.clock_analysis, f, indent=2, default=str)


def main():
    parser = argparse.ArgumentParser(description="Advanced phylogenetic analysis with temporal features")
    parser.add_argument("alignment", help="Input alignment file")
    parser.add_argument("output_tree", help="Output tree file")
    parser.add_argument("--method", choices=['iqtree', 'fasttree', 'raxml', 'neighbor_joining', 'maximum_parsimony', 'auto'],
                       default='auto', help="Tree building method")
    parser.add_argument("--model", default='auto', help="Evolutionary model")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--metadata", help="Metadata file for temporal analysis")
    parser.add_argument("--temporal-analysis", action='store_true',
                       help="Perform temporal/molecular clock analysis")
    parser.add_argument("--stats-output", help="JSON file to save tree statistics")
    parser.add_argument("--temporal-output-dir", help="Directory for temporal analysis output")
    
    args = parser.parse_args()
    
    try:
        # Build phylogenetic tree
        analyzer = AdvancedPhylogeneticAnalyzer(method=args.method, threads=args.threads)
        analyzer.build_tree(args.alignment, args.output_tree, args.metadata, args.model)
        
        # Save tree statistics
        if args.stats_output:
            combined_stats = {
                'tree_stats': analyzer.tree_stats,
                'temporal_stats': analyzer.temporal_stats,
                'method': args.method,
                'model': args.model
            }
            with open(args.stats_output, 'w') as f:
                json.dump(combined_stats, f, indent=2, default=str)
                
        # Perform temporal analysis if requested
        if args.temporal_analysis and args.metadata and args.temporal_output_dir:
            temporal_analyzer = TemporalPhylogeneticAnalyzer(threads=args.threads)
            temporal_analyzer.estimate_molecular_clock(
                args.output_tree, args.metadata, args.temporal_output_dir
            )
            
        logger.info("Phylogenetic analysis completed successfully")
        
    except Exception as e:
        logger.error(f"Phylogenetic analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
