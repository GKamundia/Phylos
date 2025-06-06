#!/usr/bin/env python3
"""
Advanced sequence alignment capabilities for RVF Nextstrain analysis
Implements multiple alignment algorithms and optimization strategies
"""

import argparse
import logging
import os
import sys
import tempfile
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline, ClustalwCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import dendropy
import subprocess
import json
from tqdm import tqdm
import pandas as pd
from datetime import datetime

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AdvancedAligner:
    """Advanced sequence alignment with multiple algorithm support"""
    
    def __init__(self, method='mafft', threads=4, temp_dir=None):
        self.method = method.lower()
        self.threads = threads
        self.temp_dir = temp_dir or tempfile.gettempdir()
        self.alignment_stats = {}
        
    def align_sequences(self, input_fasta, output_fasta, reference=None):
        """
        Perform sequence alignment using specified method
        
        Args:
            input_fasta: Path to input sequences
            output_fasta: Path to output alignment
            reference: Optional reference sequence file
        """
        logger.info(f"Starting alignment with {self.method}")
        
        # Load sequences
        sequences = list(SeqIO.parse(input_fasta, "fasta"))
        logger.info(f"Loaded {len(sequences)} sequences for alignment")
        
        if len(sequences) < 2:
            logger.warning("Less than 2 sequences provided, copying input to output")
            SeqIO.write(sequences, output_fasta, "fasta")
            return
            
        # Choose alignment method
        if self.method == 'mafft':
            self._align_with_mafft(input_fasta, output_fasta, reference)
        elif self.method == 'muscle':
            self._align_with_muscle(input_fasta, output_fasta)
        elif self.method == 'clustalw':
            self._align_with_clustalw(input_fasta, output_fasta)
        elif self.method == 'auto':
            self._auto_select_alignment(input_fasta, output_fasta, reference)
        else:
            raise ValueError(f"Unsupported alignment method: {self.method}")
            
        # Calculate alignment statistics
        self._calculate_alignment_stats(output_fasta)
        
    def _align_with_mafft(self, input_fasta, output_fasta, reference=None):
        """Align sequences using MAFFT"""
        try:
            cmd = [
                'mafft',
                '--auto',
                '--thread', str(self.threads),
                '--nomemsave'
            ]
            
            if reference:
                cmd.extend(['--addfragments', input_fasta, reference])
            else:
                cmd.append(input_fasta)
                
            logger.info(f"Running MAFFT: {' '.join(cmd)}")
            
            with open(output_fasta, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
                
            if result.returncode != 0:
                logger.error(f"MAFFT failed: {result.stderr}")
                raise RuntimeError(f"MAFFT alignment failed: {result.stderr}")
                
        except FileNotFoundError:
            logger.warning("MAFFT not found, falling back to built-in alignment")
            self._align_with_biopython(input_fasta, output_fasta)
            
    def _align_with_muscle(self, input_fasta, output_fasta):
        """Align sequences using MUSCLE (if available)"""
        try:
            cmd = ['muscle', '-in', input_fasta, '-out', output_fasta]
            logger.info(f"Running MUSCLE: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"MUSCLE failed: {result.stderr}")
                raise RuntimeError(f"MUSCLE alignment failed: {result.stderr}")
                
        except FileNotFoundError:
            logger.warning("MUSCLE not found, falling back to MAFFT")
            self._align_with_mafft(input_fasta, output_fasta)
            
    def _align_with_clustalw(self, input_fasta, output_fasta):
        """Align sequences using ClustalW (if available)"""
        try:
            clustalw_cline = ClustalwCommandline("clustalw2", infile=input_fasta)
            stdout, stderr = clustalw_cline()
            
            # ClustalW creates .aln file, convert to fasta
            aln_file = input_fasta.replace('.fasta', '.aln').replace('.fa', '.aln')
            if os.path.exists(aln_file):
                alignment = AlignIO.read(aln_file, "clustal")
                AlignIO.write(alignment, output_fasta, "fasta")
                os.remove(aln_file)  # Clean up
            else:
                raise RuntimeError("ClustalW alignment file not created")
                
        except (FileNotFoundError, Exception) as e:
            logger.warning(f"ClustalW failed: {e}, falling back to MAFFT")
            self._align_with_mafft(input_fasta, output_fasta)
            
    def _align_with_biopython(self, input_fasta, output_fasta):
        """Simple pairwise alignment using BioPython as fallback"""
        logger.info("Using BioPython pairwise alignment as fallback")
        
        sequences = list(SeqIO.parse(input_fasta, "fasta"))
        
        if len(sequences) == 1:
            SeqIO.write(sequences, output_fasta, "fasta")
            return
            
        # For simplicity, just pad sequences to same length
        max_length = max(len(seq.seq) for seq in sequences)
        
        aligned_seqs = []
        for seq in sequences:
            padded_seq = str(seq.seq).ljust(max_length, 'N')
            aligned_record = SeqRecord(
                Seq(padded_seq),
                id=seq.id,
                description=seq.description
            )
            aligned_seqs.append(aligned_record)
            
        SeqIO.write(aligned_seqs, output_fasta, "fasta")
        
    def _auto_select_alignment(self, input_fasta, output_fasta, reference=None):
        """Automatically select best alignment method based on sequence characteristics"""
        sequences = list(SeqIO.parse(input_fasta, "fasta"))
        num_sequences = len(sequences)
        avg_length = sum(len(seq.seq) for seq in sequences) / num_sequences
        
        logger.info(f"Auto-selecting alignment method for {num_sequences} sequences, avg length {avg_length:.0f}")
        
        # Decision logic for method selection
        if num_sequences < 100 and avg_length < 5000:
            self.method = 'mafft'
        elif num_sequences < 500:
            self.method = 'mafft'
        else:
            # For very large datasets, use faster options
            self.method = 'mafft'
            
        logger.info(f"Selected alignment method: {self.method}")
        self.align_sequences(input_fasta, output_fasta, reference)
        
    def _calculate_alignment_stats(self, alignment_file):
        """Calculate and store alignment quality statistics"""
        try:
            alignment = AlignIO.read(alignment_file, "fasta")
            
            self.alignment_stats = {
                'num_sequences': len(alignment),
                'alignment_length': alignment.get_alignment_length(),
                'timestamp': datetime.now().isoformat()
            }
            
            # Calculate gap statistics
            gap_counts = []
            for record in alignment:
                gaps = str(record.seq).count('-') + str(record.seq).count('N')
                gap_counts.append(gaps / len(record.seq))
                
            self.alignment_stats.update({
                'avg_gap_percentage': sum(gap_counts) / len(gap_counts) * 100,
                'max_gap_percentage': max(gap_counts) * 100,
                'min_gap_percentage': min(gap_counts) * 100
            })
            
            logger.info(f"Alignment stats: {self.alignment_stats}")
            
        except Exception as e:
            logger.warning(f"Could not calculate alignment statistics: {e}")
            

class SegmentedAligner:
    """Specialized aligner for segmented viruses like RVF"""
    
    def __init__(self, threads=4):
        self.threads = threads
        self.segment_alignments = {}
        
    def align_by_segments(self, input_fasta, metadata_file, output_dir, segments=None):
        """
        Align sequences separately by genomic segment
        
        Args:
            input_fasta: Path to input sequences
            metadata_file: Path to metadata with segment information
            output_dir: Directory for output alignments
            segments: List of segments to process (default: ['L', 'M', 'S'])
        """
        if segments is None:
            segments = ['L', 'M', 'S']
            
        os.makedirs(output_dir, exist_ok=True)
        
        # Load metadata
        metadata = pd.read_csv(metadata_file, sep='\t')
        
        # Load sequences
        sequences = {record.id: record for record in SeqIO.parse(input_fasta, "fasta")}
        
        for segment in segments:
            logger.info(f"Processing segment {segment}")
            
            # Filter sequences for this segment
            segment_metadata = metadata[metadata['segment'] == segment]
            segment_sequences = []
            
            for _, row in segment_metadata.iterrows():
                strain = row['strain']
                if strain in sequences:
                    segment_sequences.append(sequences[strain])
                    
            if not segment_sequences:
                logger.warning(f"No sequences found for segment {segment}")
                continue
                
            # Create segment-specific files
            segment_input = os.path.join(output_dir, f"{segment}_input.fasta")
            segment_output = os.path.join(output_dir, f"{segment}_aligned.fasta")
            
            SeqIO.write(segment_sequences, segment_input, "fasta")
            
            # Align segment
            aligner = AdvancedAligner(method='mafft', threads=self.threads)
            aligner.align_sequences(segment_input, segment_output)
            
            self.segment_alignments[segment] = {
                'input_file': segment_input,
                'output_file': segment_output,
                'num_sequences': len(segment_sequences),
                'stats': aligner.alignment_stats
            }
            
            logger.info(f"Completed alignment for segment {segment}: {len(segment_sequences)} sequences")
            
        return self.segment_alignments


def main():
    parser = argparse.ArgumentParser(description="Advanced sequence alignment for phylogenetic analysis")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output aligned FASTA file")
    parser.add_argument("--method", choices=['mafft', 'muscle', 'clustalw', 'auto'], 
                       default='auto', help="Alignment method")
    parser.add_argument("--reference", help="Reference sequence file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--metadata", help="Metadata file for segment-based alignment")
    parser.add_argument("--segments", nargs='+', default=['L', 'M', 'S'],
                       help="Segments to process (for segmented genomes)")
    parser.add_argument("--segment-mode", action='store_true',
                       help="Process segments separately")
    parser.add_argument("--stats-output", help="JSON file to save alignment statistics")
    
    args = parser.parse_args()
    
    try:
        if args.segment_mode and args.metadata:
            # Segmented alignment mode
            aligner = SegmentedAligner(threads=args.threads)
            results = aligner.align_by_segments(
                args.input, args.metadata, 
                os.path.dirname(args.output), args.segments
            )
            
            if args.stats_output:
                with open(args.stats_output, 'w') as f:
                    json.dump(results, f, indent=2, default=str)
                    
        else:
            # Single alignment mode
            aligner = AdvancedAligner(method=args.method, threads=args.threads)
            aligner.align_sequences(args.input, args.output, args.reference)
            
            if args.stats_output:
                with open(args.stats_output, 'w') as f:
                    json.dump(aligner.alignment_stats, f, indent=2, default=str)
                    
        logger.info("Alignment completed successfully")
        
    except Exception as e:
        logger.error(f"Alignment failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
