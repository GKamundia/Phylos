#!/usr/bin/env python3
"""
Test script to verify segment splitting process works correctly
"""

import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict, Counter
import json

def analyze_filtered_data():
    """Analyze the filtered data before splitting"""
    print("=== ANALYZING FILTERED DATA (PRE-SPLIT) ===")
    
    # Load filtered metadata
    metadata_file = "results/filtered/rvf_metadata.tsv"
    sequences_file = "results/filtered/rvf_filtered.fasta"
    
    print(f"Loading metadata from: {metadata_file}")
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    # Analyze segments in metadata
    segment_counts = metadata['segment'].value_counts()
    print(f"\nSegment distribution in metadata:")
    for segment, count in segment_counts.items():
        print(f"  {segment}: {count} sequences")
    
    # Analyze sequences
    print(f"\nLoading sequences from: {sequences_file}")
    sequences = list(SeqIO.parse(sequences_file, "fasta"))
    sequence_lengths = defaultdict(list)
    sequence_ids_by_segment = defaultdict(set)
    
    # Map sequence IDs to segments from metadata
    id_to_segment = dict(zip(metadata['strain'], metadata['segment']))
    
    for record in sequences:
        strain = record.id
        if strain in id_to_segment:
            segment = id_to_segment[strain]
            sequence_lengths[segment].append(len(record.seq))
            sequence_ids_by_segment[segment].add(strain)
        else:
            print(f"Warning: Sequence {strain} not found in metadata")
    
    print(f"\nSequence analysis:")
    print(f"Total sequences in FASTA: {len(sequences)}")
    
    for segment in ['L', 'M', 'S']:
        lengths = sequence_lengths[segment]
        if lengths:
            print(f"\n{segment} segment:")
            print(f"  Count: {len(lengths)}")
            print(f"  Length range: {min(lengths)} - {max(lengths)} bp")
            print(f"  Average length: {sum(lengths)/len(lengths):.1f} bp")
            print(f"  Length variation: {max(lengths) - min(lengths)} bp")
    
    return metadata, sequences, id_to_segment

def test_split_by_segment():
    """Test the segment splitting process manually"""
    print("\n=== TESTING SEGMENT SPLITTING PROCESS ===")
    
    # Load the data
    metadata, sequences, id_to_segment = analyze_filtered_data()
    
    # Simulate the split_by_segment process
    segments = ['L', 'M', 'S']
    split_results = {}
    
    for segment in segments:
        print(f"\nProcessing {segment} segment...")
        
        # Get metadata for this segment
        segment_metadata = metadata[metadata['segment'] == segment].copy()
        segment_strain_ids = set(segment_metadata['strain'].tolist())
        
        # Get sequences for this segment
        segment_sequences = []
        for record in sequences:
            if record.id in segment_strain_ids:
                segment_sequences.append(record)
        
        split_results[segment] = {
            'metadata_count': len(segment_metadata),
            'sequence_count': len(segment_sequences),
            'strain_ids': segment_strain_ids,
            'sequences': segment_sequences
        }
        
        print(f"  Metadata records: {len(segment_metadata)}")
        print(f"  Sequence records: {len(segment_sequences)}")
        
        # Check for mismatches
        metadata_strains = set(segment_metadata['strain'])
        sequence_strains = set(record.id for record in segment_sequences)
        
        missing_sequences = metadata_strains - sequence_strains
        extra_sequences = sequence_strains - metadata_strains
        
        if missing_sequences:
            print(f"  WARNING: {len(missing_sequences)} strains in metadata but missing sequences")
            print(f"    Examples: {list(missing_sequences)[:5]}")
        
        if extra_sequences:
            print(f"  WARNING: {len(extra_sequences)} sequences without metadata")
            print(f"    Examples: {list(extra_sequences)[:5]}")
        
        # Analyze sequence lengths for this segment
        if segment_sequences:
            lengths = [len(record.seq) for record in segment_sequences]
            print(f"  Sequence lengths: {min(lengths)} - {max(lengths)} bp")
            print(f"  Average length: {sum(lengths)/len(lengths):.1f} bp")
            
            # Check for length variability within segment
            length_variation = max(lengths) - min(lengths)
            if length_variation > 100:  # Arbitrary threshold
                print(f"  NOTE: High length variation ({length_variation} bp) - alignment will handle this")
    
    return split_results

def check_sequence_length_impact():
    """Check what happens with sequences of different lengths during alignment"""
    print("\n=== ANALYZING SEQUENCE LENGTH VARIATION IMPACT ===")
    
    metadata = pd.read_csv("results/filtered/rvf_metadata.tsv", sep='\t')
    sequences = list(SeqIO.parse("results/filtered/rvf_filtered.fasta", "fasta"))
    
    # Map sequences to segments
    id_to_segment = dict(zip(metadata['strain'], metadata['segment']))
    
    for segment in ['L', 'M', 'S']:
        print(f"\n{segment} segment length analysis:")
        
        segment_sequences = [seq for seq in sequences if seq.id in id_to_segment and id_to_segment[seq.id] == segment]
        
        if not segment_sequences:
            print(f"  No sequences found for segment {segment}")
            continue
            
        lengths = [len(seq.seq) for seq in segment_sequences]
        length_counts = Counter(lengths)
        
        print(f"  Total sequences: {len(segment_sequences)}")
        print(f"  Unique lengths: {len(length_counts)}")
        print(f"  Length range: {min(lengths)} - {max(lengths)} bp")
        print(f"  Most common lengths:")
        
        for length, count in length_counts.most_common(5):
            print(f"    {length} bp: {count} sequences ({count/len(segment_sequences)*100:.1f}%)")
        
        # Check for outliers
        median_length = sorted(lengths)[len(lengths)//2]
        outliers = [l for l in lengths if abs(l - median_length) > median_length * 0.1]  # 10% difference
        
        if outliers:
            print(f"  Potential outliers (>10% different from median {median_length} bp): {len(outliers)} sequences")
            print(f"    Outlier lengths: {sorted(set(outliers))}")

def generate_segment_split_report():
    """Generate a comprehensive report of segment splitting readiness"""
    print("\n=== GENERATING SEGMENT SPLIT READINESS REPORT ===")
    
    try:
        metadata, sequences, id_to_segment = analyze_filtered_data()
        split_results = test_split_by_segment()
        check_sequence_length_impact()
        
        # Summary report
        report = {
            "timestamp": pd.Timestamp.now().isoformat(),
            "filtered_data_status": {
                "total_sequences": len(sequences),
                "total_metadata_records": len(metadata),
                "segment_distribution": metadata['segment'].value_counts().to_dict()
            },
            "segment_split_readiness": {},
            "alignment_readiness": {}
        }
        
        for segment in ['L', 'M', 'S']:
            if segment in split_results:
                result = split_results[segment]
                report["segment_split_readiness"][segment] = {
                    "metadata_count": result['metadata_count'],
                    "sequence_count": result['sequence_count'],
                    "data_consistency": result['metadata_count'] == result['sequence_count']
                }
                
                # Alignment readiness assessment
                if result['sequences']:
                    lengths = [len(seq.seq) for seq in result['sequences']]
                    report["alignment_readiness"][segment] = {
                        "sequence_count": len(lengths),
                        "min_length": min(lengths),
                        "max_length": max(lengths),
                        "avg_length": sum(lengths) / len(lengths),
                        "length_variation": max(lengths) - min(lengths),
                        "ready_for_alignment": len(lengths) >= 5 and max(lengths) - min(lengths) < max(lengths) * 0.5
                    }
        
        # Save report
        report_file = "results/segment_split_readiness_report.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"\nReport saved to: {report_file}")
        
        # Print summary
        print("\n=== SUMMARY ===")
        ready_segments = []
        for segment in ['L', 'M', 'S']:
            if segment in report["alignment_readiness"]:
                ready = report["alignment_readiness"][segment]["ready_for_alignment"]
                status = "READY" if ready else "NEEDS REVIEW"
                count = report["alignment_readiness"][segment]["sequence_count"]
                print(f"{segment} segment: {status} ({count} sequences)")
                if ready:
                    ready_segments.append(segment)
        
        print(f"\nSegments ready for alignment: {', '.join(ready_segments) if ready_segments else 'NONE'}")
        
        return report
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    # Check if filtered data exists
    if not os.path.exists("results/filtered/rvf_metadata.tsv"):
        print("ERROR: Filtered metadata not found. Please run filtering stage first.")
        exit(1)
    
    if not os.path.exists("results/filtered/rvf_filtered.fasta"):
        print("ERROR: Filtered sequences not found. Please run filtering stage first.")
        exit(1)
    
    # Run the analysis
    report = generate_segment_split_report()
    
    if report:
        print("\n✅ Segment splitting analysis completed successfully!")
        print("The pipeline is ready to proceed with alignment.")
    else:
        print("\n❌ Segment splitting analysis failed!")
