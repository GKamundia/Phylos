#!/usr/bin/env python3
"""
Run sequence alignment using augur align in a Windows-compatible way
"""

import os
import sys
import subprocess
import shutil
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Run alignment for RVF sequences")
    parser.add_argument("--method", default="mafft", choices=["mafft", "clustal"], 
                        help="Alignment method to use")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads to use for alignment")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Define paths
    input_sequences = "results/filtered/rvf_filtered.fasta"
    output_alignment = "results/aligned/rvf_aligned.fasta"
    log_file = "logs/align_rvf.log"
    
    # Create output directories
    os.makedirs("results/aligned", exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    
    # Check if the input file exists
    if not os.path.exists(input_sequences):
        print(f"Error: Input file {input_sequences} not found")
        sys.exit(1)
    
    # Construct the augur align command
    align_cmd = [
        "augur", "align",
        "--sequences", input_sequences,
        "--output", output_alignment,
        "--nthreads", str(args.threads),
        "--method", args.method
    ]
    
    if args.method == "mafft":
        align_cmd.extend(["--mafft-options", "--nomemsave --auto"])
    
    # Run the alignment
    print(f"Running alignment with {args.method} using {args.threads} threads")
    print(f"Input: {input_sequences}")
    print(f"Output: {output_alignment}")
    
    try:
        with open(log_file, "w") as log:
            process = subprocess.run(align_cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
            
            if process.returncode == 0:
                print("Alignment completed successfully")
                print(f"Aligned sequences written to {output_alignment}")
            else:
                print(f"Alignment failed with exit code {process.returncode}")
                print(f"Check {log_file} for details")
                sys.exit(1)
    
    except Exception as e:
        print(f"Error running alignment: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
