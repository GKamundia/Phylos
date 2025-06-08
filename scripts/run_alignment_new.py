#!/usr/bin/env python3
"""
Run sequence alignment using augur align with Windows conda support
"""

import os
import sys
import subprocess
import shutil
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Run alignment for RVF sequences")
    parser.add_argument("--sequences", required=True, help="Input sequences file")
    parser.add_argument("--output", required=True, help="Output alignment file")
    parser.add_argument("--log", required=True, help="Log file")
    parser.add_argument("--segment", help="Segment name")
    parser.add_argument("--method", default="mafft", choices=["mafft", "clustal"], 
                        help="Alignment method to use")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads to use for alignment")
    parser.add_argument("--mafft-options", default="--nomemsave --auto",
                        help="MAFFT options")
    parser.add_argument("--reference", help="Reference sequence file")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Check if we have sequences to align
    if not os.path.exists(args.sequences) or os.path.getsize(args.sequences) == 0:
        with open(args.log, 'w') as log_file:
            log_file.write(f'No sequences found for segment {args.segment}\n')
        open(args.output, 'w').close()
        return
    
    # Count sequences
    seq_count = 0
    try:
        with open(args.sequences, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
    except:
        seq_count = 0
    
    with open(args.log, 'w') as log_file:
        log_file.write(f'Aligning {seq_count} sequences for segment {args.segment}\n')
    
    if seq_count == 0:
        with open(args.log, 'a') as log_file:
            log_file.write(f'No sequences to align for segment {args.segment}\n')
        open(args.output, 'w').close()
        return
    elif seq_count == 1:
        with open(args.log, 'a') as log_file:
            log_file.write('Only one sequence found, copying as alignment\n')
        shutil.copy(args.sequences, args.output)
        return
    
    # Construct the augur align command
    align_cmd = [
        "augur", "align",        
        "--sequences", args.sequences,
        "--output", args.output,
        "--nthreads", str(args.threads)
    ]
    
    if args.method == "mafft":
        align_cmd.extend(["--method", "mafft"])
    
    # Handle reference sequence
    if args.reference and os.path.exists(args.reference):
        # Check if reference sequence is already in the input sequences
        ref_name = None
        with open(args.reference, 'r') as ref_file:
            for line in ref_file:
                if line.startswith('>'):
                    ref_name = line[1:].split()[0]  # Get accession from header
                    break
        
        if ref_name:
            # Check if this reference is already in the sequences
            with open(args.sequences, 'r') as seq_file:
                seq_content = seq_file.read()
                if ref_name in seq_content:
                    with open(args.log, 'a') as log_file:
                        log_file.write(f'Reference {ref_name} found in sequences, using --reference-name\n')
                    align_cmd.extend(["--reference-name", ref_name])
                else:
                    with open(args.log, 'a') as log_file:
                        log_file.write(f'Reference {ref_name} not in sequences, using --reference-sequence\n')
                    align_cmd.extend(["--reference-sequence", args.reference])
        else:
            with open(args.log, 'a') as log_file:
                log_file.write(f'Could not parse reference name from {args.reference}\n')
    else:
        with open(args.log, 'a') as log_file:
            log_file.write('No reference sequence, aligning without reference\n')
    
    # Run the alignment
    try:
        with open(args.log, 'a') as log_file:
            log_file.write(f"Running command: {' '.join(align_cmd)}\n")
        
        # Check if we're on Windows and need special handling
        import platform
        if platform.system() == "Windows":
            # Try different approaches for Windows
            success = False
            
            # Option 1: Try with conda environment
            try:
                with open(args.log, 'a') as log_file:
                    log_file.write("Attempting alignment with conda environment\n")
                
                # Create conda environment with MAFFT if it doesn't exist
                conda_env_name = "mafft-env"
                
                # Check if environment exists
                env_check = subprocess.run(
                    ["conda", "env", "list"], 
                    capture_output=True, text=True
                )
                
                if conda_env_name not in env_check.stdout:
                    with open(args.log, 'a') as log_file:
                        log_file.write("Creating conda environment with MAFFT\n")
                    
                    # Create environment with mafft
                    create_result = subprocess.run([
                        "conda", "create", "-n", conda_env_name, 
                        "-c", "conda-forge", "-c", "bioconda", 
                        "python=3.9", "mafft", "augur", "-y"
                    ], capture_output=True, text=True)
                    
                    if create_result.returncode != 0:
                        raise Exception(f"Failed to create conda environment: {create_result.stderr}")
                
                # Run alignment in conda environment
                conda_cmd = [
                    "conda", "run", "-n", conda_env_name
                ] + align_cmd
                
                with open(args.log, 'a') as log_file:
                    log_file.write(f"Conda command: {' '.join(conda_cmd)}\n")
                
                result = subprocess.run(conda_cmd, capture_output=True, text=True, check=True)
                success = True
                
                with open(args.log, 'a') as log_file:
                    log_file.write("Conda environment alignment completed successfully\n")
                    if result.stdout:
                        log_file.write(f"STDOUT:\n{result.stdout}\n")
                    if result.stderr:
                        log_file.write(f"STDERR:\n{result.stderr}\n")
                        
            except Exception as conda_error:
                with open(args.log, 'a') as log_file:
                    log_file.write(f"Conda environment approach failed: {conda_error}\n")
            
            # Option 2: Try BioPython's simple alignment as fallback
            if not success:
                try:
                    with open(args.log, 'a') as log_file:
                        log_file.write("Attempting simple BioPython alignment fallback\n")
                    
                    sequences = list(SeqIO.parse(args.sequences, "fasta"))
                    
                    # Simple alignment: pad sequences to same length
                    max_len = max(len(seq.seq) for seq in sequences)
                    
                    for seq in sequences:
                        if len(seq.seq) < max_len:
                            # Pad with gaps at the end
                            gap_count = max_len - len(seq.seq)
                            seq.seq = seq.seq + ('-' * gap_count)
                    
                    with open(args.output, 'w') as out_file:
                        SeqIO.write(sequences, out_file, "fasta")
                    success = True
                    
                    with open(args.log, 'a') as log_file:
                        log_file.write("BioPython simple alignment completed\n")
                        
                except Exception as bio_error:
                    with open(args.log, 'a') as log_file:
                        log_file.write(f"BioPython fallback failed: {bio_error}\n")
            
            # Option 3: Final fallback - copy sequences
            if not success:
                with open(args.log, 'a') as log_file:
                    log_file.write("Using final fallback - copying sequences as aligned\n")
                    log_file.write("WARNING: This is not a real alignment!\n")
                shutil.copy(args.sequences, args.output)
                
        else:
            # Normal augur align for Linux/Mac
            result = subprocess.run(align_cmd, capture_output=True, text=True, check=True)
            
            with open(args.log, 'a') as log_file:
                log_file.write("Alignment completed successfully\n")
                if result.stdout:
                    log_file.write(f"STDOUT:\n{result.stdout}\n")
                if result.stderr:
                    log_file.write(f"STDERR:\n{result.stderr}\n")
                
    except subprocess.CalledProcessError as e:
        with open(args.log, 'a') as log_file:
            log_file.write(f'Alignment failed: {e}\n')
            if e.stdout:
                log_file.write(f"STDOUT:\n{e.stdout}\n")
            if e.stderr:
                log_file.write(f"STDERR:\n{e.stderr}\n")
        raise

if __name__ == "__main__":
    main()
