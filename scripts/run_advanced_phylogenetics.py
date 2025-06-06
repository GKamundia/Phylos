#!/usr/bin/env python3
"""
Integration script for advanced phylogenetic analysis pipeline
Demonstrates how to run the complete workflow with the new advanced features
"""

import argparse
import logging
import os
import sys
import subprocess
from pathlib import Path
import json

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_advanced_phylogenetic_pipeline(sequences_file, metadata_file, output_dir, 
                                      alignment_method='auto', tree_method='auto', 
                                      threads=4, include_temporal=True):
    """
    Run the complete advanced phylogenetic analysis pipeline
    
    Args:
        sequences_file: Path to input FASTA sequences
        metadata_file: Path to metadata TSV file
        output_dir: Output directory for all results
        alignment_method: Alignment method ('auto', 'mafft', 'muscle', 'clustalw')
        tree_method: Tree building method ('auto', 'iqtree', 'fasttree', etc.)
        threads: Number of threads to use
        include_temporal: Whether to perform temporal analysis
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output files
    aligned_file = os.path.join(output_dir, 'aligned_sequences.fasta')
    alignment_stats = os.path.join(output_dir, 'alignment_stats.json')
    tree_file = os.path.join(output_dir, 'phylogenetic_tree.nwk')
    tree_stats = os.path.join(output_dir, 'tree_stats.json')
    temporal_dir = os.path.join(output_dir, 'temporal_analysis')
    
    pipeline_results = {
        'input_files': {
            'sequences': sequences_file,
            'metadata': metadata_file
        },
        'parameters': {
            'alignment_method': alignment_method,
            'tree_method': tree_method,
            'threads': threads,
            'include_temporal': include_temporal
        },
        'output_files': {},
        'steps_completed': [],
        'errors': []
    }
    
    try:
        # Step 1: Advanced Alignment
        logger.info("Step 1: Running advanced sequence alignment")
        
        alignment_cmd = [
            'python', 'scripts/advanced_alignment.py',
            sequences_file,
            aligned_file,
            '--method', alignment_method,
            '--threads', str(threads),
            '--stats-output', alignment_stats
        ]
        
        result = subprocess.run(alignment_cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            pipeline_results['steps_completed'].append('alignment')
            pipeline_results['output_files']['alignment'] = aligned_file
            pipeline_results['output_files']['alignment_stats'] = alignment_stats
            logger.info("✓ Alignment completed successfully")
        else:
            raise RuntimeError(f"Alignment failed: {result.stderr}")
            
        # Step 2: Advanced Tree Building
        logger.info("Step 2: Building phylogenetic tree")
        
        tree_cmd = [
            'python', 'scripts/advanced_phylogenetics.py',
            aligned_file,
            tree_file,
            '--method', tree_method,
            '--threads', str(threads),
            '--metadata', metadata_file,
            '--stats-output', tree_stats
        ]
        
        if include_temporal:
            tree_cmd.extend(['--temporal-analysis', '--temporal-output-dir', temporal_dir])
        
        result = subprocess.run(tree_cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            pipeline_results['steps_completed'].append('tree_building')
            pipeline_results['output_files']['tree'] = tree_file
            pipeline_results['output_files']['tree_stats'] = tree_stats
            logger.info("✓ Tree building completed successfully")
        else:
            raise RuntimeError(f"Tree building failed: {result.stderr}")
            
        # Step 3: Temporal Analysis (if requested)
        if include_temporal:
            logger.info("Step 3: Running temporal analysis")
            
            temporal_cmd = [
                'python', 'scripts/temporal_analysis.py',
                metadata_file,
                '--tree', tree_file,
                '--output-dir', temporal_dir
            ]
            
            result = subprocess.run(temporal_cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                pipeline_results['steps_completed'].append('temporal_analysis')
                pipeline_results['output_files']['temporal_dir'] = temporal_dir
                pipeline_results['output_files']['temporal_summary'] = os.path.join(temporal_dir, 'temporal_analysis_summary.txt')
                logger.info("✓ Temporal analysis completed successfully")
            else:
                logger.warning(f"Temporal analysis failed: {result.stderr}")
                pipeline_results['errors'].append(f"Temporal analysis failed: {result.stderr}")
        
        # Step 4: Generate Summary Report
        logger.info("Step 4: Generating summary report")
        
        summary_file = os.path.join(output_dir, 'pipeline_summary.json')
        with open(summary_file, 'w') as f:
            json.dump(pipeline_results, f, indent=2, default=str)
            
        pipeline_results['output_files']['summary'] = summary_file
        
        # Create human-readable summary
        summary_txt = os.path.join(output_dir, 'pipeline_summary.txt')
        with open(summary_txt, 'w') as f:
            f.write("ADVANCED PHYLOGENETIC ANALYSIS PIPELINE SUMMARY\n")
            f.write("=" * 55 + "\n\n")
            
            f.write("INPUT FILES:\n")
            f.write(f"  Sequences: {sequences_file}\n")
            f.write(f"  Metadata: {metadata_file}\n\n")
            
            f.write("PARAMETERS:\n")
            f.write(f"  Alignment method: {alignment_method}\n")
            f.write(f"  Tree method: {tree_method}\n")
            f.write(f"  Threads: {threads}\n")
            f.write(f"  Temporal analysis: {include_temporal}\n\n")
            
            f.write("COMPLETED STEPS:\n")
            for step in pipeline_results['steps_completed']:
                f.write(f"  ✓ {step.replace('_', ' ').title()}\n")
            f.write("\n")
            
            if pipeline_results['errors']:
                f.write("ERRORS/WARNINGS:\n")
                for error in pipeline_results['errors']:
                    f.write(f"  ⚠ {error}\n")
                f.write("\n")
            
            f.write("OUTPUT FILES:\n")
            for file_type, file_path in pipeline_results['output_files'].items():
                if os.path.exists(file_path):
                    f.write(f"  {file_type}: {file_path}\n")
                    
        pipeline_results['output_files']['summary_txt'] = summary_txt
        
        logger.info("✓ Pipeline completed successfully")
        logger.info(f"Results saved to: {output_dir}")
        
        return pipeline_results
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        pipeline_results['errors'].append(str(e))
        pipeline_results['status'] = 'failed'
        
        # Save error report
        error_file = os.path.join(output_dir, 'pipeline_error.json')
        with open(error_file, 'w') as f:
            json.dump(pipeline_results, f, indent=2, default=str)
            
        return pipeline_results


def main():
    parser = argparse.ArgumentParser(description="Advanced phylogenetic analysis pipeline")
    parser.add_argument("sequences", help="Input FASTA sequences file")
    parser.add_argument("metadata", help="Input metadata TSV file") 
    parser.add_argument("output_dir", help="Output directory for results")
    parser.add_argument("--alignment-method", choices=['auto', 'mafft', 'muscle', 'clustalw'],
                       default='auto', help="Alignment method")
    parser.add_argument("--tree-method", choices=['auto', 'iqtree', 'fasttree', 'raxml', 
                                                 'neighbor_joining', 'maximum_parsimony'],
                       default='auto', help="Tree building method")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--no-temporal", action='store_true',
                       help="Skip temporal analysis")
    parser.add_argument("--quick", action='store_true',
                       help="Use faster methods for large datasets")
    
    args = parser.parse_args()
    
    # Adjust methods for quick mode
    if args.quick:
        alignment_method = 'mafft' if args.alignment_method == 'auto' else args.alignment_method
        tree_method = 'fasttree' if args.tree_method == 'auto' else args.tree_method
    else:
        alignment_method = args.alignment_method
        tree_method = args.tree_method
    
    # Validate input files
    if not os.path.exists(args.sequences):
        logger.error(f"Sequences file not found: {args.sequences}")
        sys.exit(1)
        
    if not os.path.exists(args.metadata):
        logger.error(f"Metadata file not found: {args.metadata}")
        sys.exit(1)
    
    # Run pipeline
    results = run_advanced_phylogenetic_pipeline(
        sequences_file=args.sequences,
        metadata_file=args.metadata,
        output_dir=args.output_dir,
        alignment_method=alignment_method,
        tree_method=tree_method,
        threads=args.threads,
        include_temporal=not args.no_temporal
    )
    
    # Exit with appropriate code
    if results.get('status') == 'failed' or results.get('errors'):
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
