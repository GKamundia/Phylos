#!/usr/bin/env python3
"""
Run tree building using augur tree with Windows WSL support
"""

import os
import sys
import subprocess
import argparse
import platform

def parse_args():
    parser = argparse.ArgumentParser(description="Run tree building for sequences")
    parser.add_argument("--alignment", required=True, help="Input alignment file")
    parser.add_argument("--output", required=True, help="Output tree file")
    parser.add_argument("--log", required=True, help="Log file")
    parser.add_argument("--method", default="iqtree", choices=["iqtree", "fasttree"], 
                        help="Tree building method")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads to use")
    parser.add_argument("--substitution-model", default="GTR",
                        help="Substitution model for tree building")
    parser.add_argument("--iqtree-args", default="-ninit 2 -n 2",
                        help="Additional IQ-TREE arguments")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Check if we have alignment to build tree from
    if not os.path.exists(args.alignment) or os.path.getsize(args.alignment) == 0:
        with open(args.log, 'w') as log_file:
            log_file.write('No alignment found for tree building\n')
        open(args.output, 'w').close()
        return
    
    # Count sequences in alignment
    seq_count = 0
    try:
        with open(args.alignment, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
    except:
        seq_count = 0
    
    with open(args.log, 'w') as log_file:
        log_file.write(f'Building tree from {seq_count} aligned sequences\n')
    
    if seq_count < 2:
        with open(args.log, 'a') as log_file:
            log_file.write('Need at least 2 sequences to build a tree\n')
        open(args.output, 'w').close()
        return
    
    # Construct the augur tree command
    tree_cmd = [
        "augur", "tree",
        "--alignment", args.alignment,
        "--output", args.output,
        "--nthreads", str(args.threads)
    ]
    
    if args.method == "iqtree":
        tree_cmd.extend([
            "--method", "iqtree",
            "--substitution-model", args.substitution_model
        ])
        # Add IQ-TREE specific arguments
        if args.iqtree_args:
            tree_cmd.extend(args.iqtree_args.split())
    else:
        tree_cmd.extend(["--method", "fasttree"])
    
    # Run the tree building
    try:
        with open(args.log, 'a') as log_file:
            log_file.write(f"Running command: {' '.join(tree_cmd)}\n")
        
        # Check if we're on Windows and need special handling
        if platform.system() == "Windows":
            # Try different approaches for Windows
            success = False
            
            # Option 1: Try with WSL (Windows Subsystem for Linux)
            try:
                with open(args.log, 'a') as log_file:
                    log_file.write("Attempting tree building with WSL\n")
                
                # Convert Windows paths to WSL paths
                wsl_alignment = args.alignment.replace('C:', '/mnt/c').replace('\\', '/')
                wsl_output = args.output.replace('C:', '/mnt/c').replace('\\', '/')
                
                # Check if WSL has augur in virtual environment
                check_augur = subprocess.run([
                    "wsl", "bash", "-c", ". /home/anarchy/augur-env/bin/activate && which augur"
                ], capture_output=True, text=True)
                
                if check_augur.returncode != 0:
                    raise Exception("Augur not found in WSL virtual environment")
                
                # Construct WSL command with virtual environment activation
                if args.method == "iqtree":
                    augur_cmd = f". /home/anarchy/augur-env/bin/activate && augur tree --alignment '{wsl_alignment}' --output '{wsl_output}' --method iqtree --substitution-model {args.substitution_model} --nthreads {args.threads}"
                    if args.iqtree_args:
                        augur_cmd += f" {args.iqtree_args}"
                else:
                    augur_cmd = f". /home/anarchy/augur-env/bin/activate && augur tree --alignment '{wsl_alignment}' --output '{wsl_output}' --method fasttree --nthreads {args.threads}"
                
                wsl_tree_cmd = ["wsl", "bash", "-c", augur_cmd]
                
                with open(args.log, 'a') as log_file:
                    log_file.write(f"WSL command: {' '.join(wsl_tree_cmd)}\n")
                
                result = subprocess.run(wsl_tree_cmd, capture_output=True, text=True, check=True)
                success = True
                
                with open(args.log, 'a') as log_file:
                    log_file.write("WSL tree building completed successfully\n")
                    if result.stdout:
                        log_file.write(f"STDOUT:\n{result.stdout}\n")
                    if result.stderr:
                        log_file.write(f"STDERR:\n{result.stderr}\n")
                        
            except Exception as wsl_error:
                with open(args.log, 'a') as log_file:
                    log_file.write(f"WSL approach failed: {wsl_error}\n")
            
            # Option 2: Try with conda environment (if WSL failed)
            if not success:
                try:
                    with open(args.log, 'a') as log_file:
                        log_file.write("Attempting tree building with conda environment\n")
                    
                    # Try to run in conda environment
                    conda_env_name = "mafft-env"
                    conda_cmd = [
                        "conda", "run", "-n", conda_env_name
                    ] + tree_cmd
                    
                    with open(args.log, 'a') as log_file:
                        log_file.write(f"Conda command: {' '.join(conda_cmd)}\n")
                    
                    result = subprocess.run(conda_cmd, capture_output=True, text=True, check=True)
                    success = True
                    
                    with open(args.log, 'a') as log_file:
                        log_file.write("Conda environment tree building completed successfully\n")
                        if result.stdout:
                            log_file.write(f"STDOUT:\n{result.stdout}\n")
                        if result.stderr:
                            log_file.write(f"STDERR:\n{result.stderr}\n")
                            
                except Exception as conda_error:
                    with open(args.log, 'a') as log_file:
                        log_file.write(f"Conda environment approach failed: {conda_error}\n")
            
            # Option 3: Final fallback - create empty tree
            if not success:
                with open(args.log, 'a') as log_file:
                    log_file.write("Using final fallback - creating empty tree file\n")
                    log_file.write("WARNING: No proper tree was built!\n")
                open(args.output, 'w').close()
                
        else:
            # Normal augur tree for Linux/Mac
            result = subprocess.run(tree_cmd, capture_output=True, text=True, check=True)
            
            with open(args.log, 'a') as log_file:
                log_file.write("Tree building completed successfully\n")
                if result.stdout:
                    log_file.write(f"STDOUT:\n{result.stdout}\n")
                if result.stderr:
                    log_file.write(f"STDERR:\n{result.stderr}\n")
                
    except subprocess.CalledProcessError as e:
        with open(args.log, 'a') as log_file:
            log_file.write(f'Tree building failed: {e}\n')
            if e.stdout:
                log_file.write(f"STDOUT:\n{e.stdout}\n")
            if e.stderr:
                log_file.write(f"STDERR:\n{e.stderr}\n")
        raise

if __name__ == "__main__":
    main()
