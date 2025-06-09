#!/usr/bin/env python3
"""
Analyze exact sequence lengths for RVF segments.

This script generates detailed length distribution analysis with exact counts
for each unique sequence length in the L, M, and S segments.

Usage:
    python scripts/analyze_exact_sequence_lengths.py

Output:
    - Console output with detailed length distributions
    - Categorized sequences (partial, near-complete, complete)
    - Quality assessment metrics
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
from collections import Counter
from datetime import datetime

def analyze_segment_lengths(segment_name, fasta_path, reference_length):
    """Analyze exact sequence lengths for a segment."""
    if not os.path.exists(fasta_path):
        print(f"❌ File not found: {fasta_path}")
        return None
    
    # Read sequences and get lengths
    lengths = []
    with open(fasta_path, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            lengths.append(len(record.seq))
    
    if not lengths:
        print(f"❌ No sequences found in {fasta_path}")
        return None
    
    # Count sequences by exact length
    length_counts = Counter(lengths)
    total_sequences = len(lengths)
    unique_lengths = len(length_counts)
    
    # Categorize sequences
    partial_count = 0
    near_complete_count = 0
    complete_count = 0
    reference_count = 0
    
    for length, count in length_counts.items():
        if length == reference_length:
            reference_count = count
            complete_count += count
        elif length >= reference_length * 0.9:  # Within 90% of reference
            if length > reference_length:
                complete_count += count
            else:
                near_complete_count += count
        else:
            partial_count += count
    
    # Generate report
    print(f"\n{'='*50}")
    print(f"{segment_name} SEGMENT EXACT LENGTH ANALYSIS")
    print(f"{'='*50}")
    print(f"Total sequences: {total_sequences}")
    print(f"Unique lengths: {unique_lengths}")
    print(f"Length range: {min(lengths):,} - {max(lengths):,} bp")
    print(f"Reference length: {reference_length:,} bp")
    
    print(f"\nQUALITY DISTRIBUTION:")
    print(f"  Complete sequences: {complete_count:>3} ({complete_count/total_sequences*100:>5.1f}%)")
    print(f"  Near-complete:      {near_complete_count:>3} ({near_complete_count/total_sequences*100:>5.1f}%)")
    print(f"  Partial sequences:  {partial_count:>3} ({partial_count/total_sequences*100:>5.1f}%)")
    print(f"  At reference length: {reference_count:>3} ({reference_count/total_sequences*100:>5.1f}%)")
    
    print(f"\nEXACT LENGTH DISTRIBUTION:")
    print(f"{'Length (bp)':>12} | {'Count':>5} | {'Percentage':>10} | {'Category':>15}")
    print(f"{'-'*12} | {'-'*5} | {'-'*10} | {'-'*15}")
    
    for length in sorted(length_counts.keys()):
        count = length_counts[length]
        percentage = (count / total_sequences) * 100
        
        # Determine category
        if length == reference_length:
            category = "Reference"
        elif length >= reference_length * 0.9:
            if length > reference_length:
                category = "Complete"
            else:
                category = "Near-complete"
        else:
            category = "Partial"
        
        marker = "**" if length == reference_length or count >= total_sequences * 0.05 else "  "
        print(f"{marker}{length:>10,} | {count:>5} | {percentage:>9.1f}% | {category:>15}")
    
    return {
        'segment': segment_name,
        'total': total_sequences,
        'unique_lengths': unique_lengths,
        'min_length': min(lengths),
        'max_length': max(lengths),
        'reference_length': reference_length,
        'complete': complete_count,
        'near_complete': near_complete_count,
        'partial': partial_count,
        'reference_exact': reference_count,
        'length_counts': length_counts
    }

def main():
    """Main analysis function."""
    print("RVF NEXTSTRAIN PIPELINE - EXACT SEQUENCE LENGTH ANALYSIS")
    print(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*60)
    
    # Define segments and their reference lengths
    segments = {
        'L': {
            'path': 'results/segments/L/raw/rvf_L_sequences.fasta',
            'reference_length': 6404
        },
        'M': {
            'path': 'results/segments/M/raw/rvf_M_sequences.fasta',
            'reference_length': 3885
        },
        'S': {
            'path': 'results/segments/S/raw/rvf_S_sequences.fasta',
            'reference_length': 1690
        }
    }
    
    results = {}
    total_sequences = 0
    total_complete = 0
    
    # Analyze each segment
    for segment_name, config in segments.items():
        result = analyze_segment_lengths(
            segment_name, 
            config['path'], 
            config['reference_length']
        )
        if result:
            results[segment_name] = result
            total_sequences += result['total']
            total_complete += result['complete'] + result['near_complete']
    
    # Overall summary
    if results:
        print(f"\n{'='*60}")
        print("OVERALL QUALITY SUMMARY")
        print(f"{'='*60}")
        print(f"Total sequences across all segments: {total_sequences:,}")
        print(f"Complete + near-complete sequences: {total_complete:,} ({total_complete/total_sequences*100:.1f}%)")
        print(f"Partial sequences: {total_sequences - total_complete:,} ({(total_sequences - total_complete)/total_sequences*100:.1f}%)")
        
        print(f"\nSEGMENT-BY-SEGMENT SUMMARY:")
        print(f"{'Segment':>8} | {'Total':>6} | {'Complete':>9} | {'Near-Complete':>13} | {'Partial':>8} | {'Quality':>8}")
        print(f"{'-'*8} | {'-'*6} | {'-'*9} | {'-'*13} | {'-'*8} | {'-'*8}")
        
        for segment_name, result in results.items():
            quality_pct = (result['complete'] + result['near_complete']) / result['total'] * 100
            print(f"{segment_name:>8} | {result['total']:>6} | {result['complete']:>9} | "
                  f"{result['near_complete']:>13} | {result['partial']:>8} | {quality_pct:>7.1f}%")
        
        print(f"\n✅ Analysis complete. Data quality: {'EXCELLENT' if total_complete/total_sequences > 0.95 else 'GOOD' if total_complete/total_sequences > 0.85 else 'ACCEPTABLE'}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
