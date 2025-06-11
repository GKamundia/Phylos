#!/usr/bin/env python3
from Bio import SeqIO
import os

# Check segment files and sequence lengths
segments = ['L', 'M', 'S']
ref_lengths = {'L': 6404, 'M': 3885, 'S': 1690}

print("=== SEGMENT SEQUENCE LENGTH VERIFICATION ===")

for seg in segments:
    seq_file = f"results/segments/{seg}/raw/rvf_{seg}_sequences.fasta"
    if os.path.exists(seq_file):
        sequences = list(SeqIO.parse(seq_file, "fasta"))
        if sequences:
            lengths = [len(record.seq) for record in sequences]
            expected = ref_lengths[seg]
            
            print(f"{seg} segment: {len(sequences)} sequences")
            print(f"  Expected length: {expected}bp")
            print(f"  Actual lengths: {set(lengths)}")
            print(f"  All correct: {'YES' if all(l == expected for l in lengths) else 'NO'}")
            
            if len(set(lengths)) > 1:
                print(f"  WARNING: Multiple lengths found!")
        else:
            print(f"{seg} segment: No sequences found")
    else:
        print(f"{seg} segment: File not found")

print("\n=== VERIFICATION COMPLETE ===")
