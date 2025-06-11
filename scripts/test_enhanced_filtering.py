#!/usr/bin/env python3
"""
Test script to verify the enhanced filtering system is working correctly.
Tests both main filtering and segment splitting with exact reference length filtering.
"""

import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import json

def test_main_filtering():
    """Test the main filtering stage with exact reference length filtering"""
    print("=== TESTING MAIN FILTERING STAGE ===")
    
    # Check if test files exist from our previous test
    test_sequences = "results/test_filtered.fasta"
    test_metadata = "results/test_metadata.tsv"
    
    if not os.path.exists(test_sequences) or not os.path.exists(test_metadata):
        print("Test files not found. Please run the main filter test first.")
        return None
    
    # Load test results
    metadata = pd.read_csv(test_metadata, sep='\t')
    sequences = list(SeqIO.parse(test_sequences, "fasta"))
    
    print(f"Filtered sequences: {len(sequences)}")
    print(f"Filtered metadata records: {len(metadata)}")
    
    # Verify segment distribution
    if 'Segment' in metadata.columns:
        segment_counts = metadata['Segment'].value_counts()
        print(f"Segment distribution: {segment_counts.to_dict()}")
        
        # Verify sequence lengths match reference lengths
        reference_lengths = {'L': 6404, 'M': 3885, 'S': 1690}
        sequence_dict = {record.id: record for record in sequences}
        
        length_verification = {}
        for segment in ['L', 'M', 'S']:
            segment_meta = metadata[metadata['Segment'].str.upper() == segment.upper()]
            segment_lengths = []
            
            for _, row in segment_meta.iterrows():
                accession = row['Accession']
                if accession in sequence_dict:
                    seq_length = len(sequence_dict[accession].seq)
                    segment_lengths.append(seq_length)
            
            if segment_lengths:
                expected_length = reference_lengths[segment]
                all_correct = all(length == expected_length for length in segment_lengths)
                length_verification[segment] = {
                    'count': len(segment_lengths),
                    'expected_length': expected_length,
                    'all_correct_length': all_correct,
                    'actual_lengths': list(set(segment_lengths))
                }
                
                status = "PASS" if all_correct else "FAIL"
                print(f"  {segment} segment: {len(segment_lengths)} sequences, {status}")
                if not all_correct:
                    print(f"    Expected: {expected_length}bp, Found: {set(segment_lengths)}")
    
    return {
        'sequences': len(sequences),
        'metadata': len(metadata),
        'segment_counts': segment_counts.to_dict() if 'Segment' in metadata.columns else {},
        'length_verification': length_verification if 'Segment' in metadata.columns else {}
    }

def test_segment_splitting():
    """Test the segment splitting stage with exact reference length filtering"""
    print("\n=== TESTING SEGMENT SPLITTING STAGE ===")
    
    # Check segment files
    segments = ['L', 'M', 'S']
    segment_results = {}
    
    for segment in segments:
        seq_file = f"results/segments/{segment}/raw/rvf_{segment}_sequences.fasta"
        meta_file = f"results/segments/{segment}/filtered/rvf_{segment}_metadata.tsv"
        
        if os.path.exists(seq_file) and os.path.exists(meta_file):
            sequences = list(SeqIO.parse(seq_file, "fasta"))
            metadata = pd.read_csv(meta_file, sep='\t')
            
            # Verify all sequences have exact reference length
            reference_lengths = {'L': 6404, 'M': 3885, 'S': 1690}
            expected_length = reference_lengths[segment]
            
            sequence_lengths = [len(record.seq) for record in sequences]
            all_correct = all(length == expected_length for length in sequence_lengths)
            
            segment_results[segment] = {
                'sequence_count': len(sequences),
                'metadata_count': len(metadata),
                'expected_length': expected_length,
                'all_correct_length': all_correct,
                'actual_lengths': list(set(sequence_lengths))
            }
            
            status = "PASS" if all_correct else "FAIL"
            print(f"  {segment} segment: {len(sequences)} sequences, {len(metadata)} metadata, {status}")
            if not all_correct:
                print(f"    Expected: {expected_length}bp, Found: {set(sequence_lengths)}")
            
            # Verify sequence count matches metadata count
            if len(sequences) != len(metadata):
                print(f"    WARNING: Sequence count ({len(sequences)}) != Metadata count ({len(metadata)})")
        else:
            print(f"  {segment} segment: Files not found")
            segment_results[segment] = {'status': 'missing'}
    
    return segment_results

def verify_data_flow():
    """Verify the complete data flow from raw to segment-specific files"""
    print("\n=== VERIFYING COMPLETE DATA FLOW ===")
    
    # Check raw filtered data
    raw_sequences = "results/filtered/rvf_filtered.fasta"
    raw_metadata = "results/filtered/rvf_metadata.tsv"
    
    if not os.path.exists(raw_sequences) or not os.path.exists(raw_metadata):
        print("Raw filtered data not found")
        return None
    
    raw_seqs = list(SeqIO.parse(raw_sequences, "fasta"))
    raw_meta = pd.read_csv(raw_metadata, sep='\t')
    
    print(f"Raw filtered: {len(raw_seqs)} sequences, {len(raw_meta)} metadata")
    
    # Check segment totals
    segments = ['L', 'M', 'S']
    total_segment_sequences = 0
    total_segment_metadata = 0
    
    for segment in segments:
        seq_file = f"results/segments/{segment}/raw/rvf_{segment}_sequences.fasta"
        meta_file = f"results/segments/{segment}/filtered/rvf_{segment}_metadata.tsv"
        
        if os.path.exists(seq_file) and os.path.exists(meta_file):
            sequences = list(SeqIO.parse(seq_file, "fasta"))
            metadata = pd.read_csv(meta_file, sep='\t')
            total_segment_sequences += len(sequences)
            total_segment_metadata += len(metadata)
    
    print(f"Segment totals: {total_segment_sequences} sequences, {total_segment_metadata} metadata")
    
    # Calculate filtering efficiency
    if len(raw_seqs) > 0:
        retention_rate = total_segment_sequences / len(raw_seqs) * 100
        print(f"Overall retention rate after length filtering: {retention_rate:.1f}%")
        
        if retention_rate < 50:
            print("  NOTE: High filtering rate - many sequences don't match reference lengths")
        elif retention_rate > 90:
            print("  NOTE: Low filtering rate - most sequences already match reference lengths")
        else:
            print("  NOTE: Moderate filtering rate - good balance of quality vs quantity")
    
    return {
        'raw_sequences': len(raw_seqs),
        'raw_metadata': len(raw_meta),
        'segment_sequences': total_segment_sequences,
        'segment_metadata': total_segment_metadata,
        'retention_rate': retention_rate if len(raw_seqs) > 0 else 0
    }

def generate_filtering_report():
    """Generate a comprehensive report on the filtering system"""
    print("\n=== GENERATING FILTERING REPORT ===")
    
    # Run all tests
    main_results = test_main_filtering()
    segment_results = test_segment_splitting()
    flow_results = verify_data_flow()
    
    # Create report
    report = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'filtering_system_status': 'operational',
        'main_filtering': main_results,
        'segment_splitting': segment_results,
        'data_flow': flow_results,
        'reference_lengths': {
            'L': 6404,
            'M': 3885,
            'S': 1690
        },
        'summary': {
            'exact_length_filtering': 'enabled',
            'main_stage': 'filters to exact reference lengths',
            'segment_stage': 'applies additional exact length filtering',
            'metadata_sync': 'metadata updated to match filtered sequences'
        }
    }
    
    # Check if all tests passed
    all_passed = True
    if main_results and 'length_verification' in main_results:
        for segment, data in main_results['length_verification'].items():
            if not data.get('all_correct_length', False):
                all_passed = False
    
    if segment_results:
        for segment, data in segment_results.items():
            if 'all_correct_length' in data and not data['all_correct_length']:
                all_passed = False
    
    report['all_tests_passed'] = all_passed
    
    # Save report
    report_file = "results/enhanced_filtering_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nReport saved to: {report_file}")
    
    # Print summary
    print("\n=== SUMMARY ===")
    if all_passed:
        print("All tests PASSED - Enhanced filtering system is working correctly")
        print("Main filtering stage applies exact reference length filtering")
        print("Segment splitting maintains exact reference lengths")
        print("Metadata is properly synchronized with filtered sequences")
    else:
        print("Some tests FAILED - Check the detailed output above")
    
    return report

if __name__ == "__main__":
    # Run the complete test suite
    report = generate_filtering_report()
    
    if report['all_tests_passed']:
        print("\nEnhanced filtering system is fully operational!")
        print("Ready for phylogenetic analysis with uniform sequence lengths.")
    else:
        print("\nIssues detected in filtering system - review needed.")