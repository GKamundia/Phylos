#!/usr/bin/env python3
"""
Run two-step IQ-TREE2 tree building with Windows WSL support and publication-ready workflow
"""

import os
import sys
import subprocess
import argparse
import platform
import re


def parse_args():
    parser = argparse.ArgumentParser(description="Run two-step IQ-TREE2 tree building for sequences")
    parser.add_argument("--alignment", required=True, help="Input alignment file")
    parser.add_argument("--output", required=True, help="Output tree file (final support tree)")
    parser.add_argument("--log", required=True, help="Log file (will contain both steps)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for reproducibility")
    parser.add_argument("--ml-prefix", default=None, help="Prefix for ML step outputs (default: output file basename)")
    parser.add_argument("--support-prefix", default=None, help="Prefix for support step outputs (default: output file basename + _support)")
    parser.add_argument("--n", type=int, default=10, help="Number of independent ML searches (Step 1)")
    parser.add_argument("--ninit", type=int, default=100, help="Number of initial parsimony trees (Step 1)")
    parser.add_argument("--bb", type=int, default=1000, help="Ultrafast bootstrap replicates (Step 2)")
    parser.add_argument("--alrt", type=int, default=1000, help="SH-aLRT support replicates (Step 2)")
    parser.add_argument("--nstop", type=int, default=100, help="Number of unsuccessful iterations to stop (Step 2, default: 100)")
    return parser.parse_args()

def run_cmd(cmd, log_file, use_wsl=False):
    # Log the command as a single line for reproducibility
    cmd_str = ' '.join(cmd)
    if use_wsl:
        wsl_cmd = ["wsl", "bash", "-c", cmd_str]
        with open(log_file, 'a') as log:
            log.write(f"\nRunning WSL command (one line): {' '.join(wsl_cmd)}\n")
        result = subprocess.run(wsl_cmd, capture_output=True, text=True)
    else:
        with open(log_file, 'a') as log:
            log.write(f"\nRunning command (one line): {cmd_str}\n")
        result = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_file, 'a') as log:
        log.write(f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}\n")
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result

def extract_best_fit_model(iqtree_file):
    # Look for 'Best-fit model according to BIC:' first, then fallback to 'Best-fit model:'
    with open(iqtree_file, 'r') as f:
        for line in f:
            if 'Best-fit model according to BIC:' in line:
                return line.strip().split(':', 1)[1].strip()
        f.seek(0)
        for line in f:
            if 'Best-fit model:' in line:
                return line.strip().split(':', 1)[1].strip()
    raise RuntimeError(f"Best-fit model not found in {iqtree_file}")

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    log_file = args.log
    # Detect if we are on Windows and need to use WSL
    use_wsl = False
    if platform.system() == "Windows":
        use_wsl = True
        # Convert paths to WSL format
        def to_wsl_path(path):
            return path.replace('C:', '/mnt/c').replace('c:', '/mnt/c').replace('\\', '/').replace('\\', '/')
        args.alignment = to_wsl_path(args.alignment)
        args.output = to_wsl_path(args.output)
        if args.ml_prefix:
            args.ml_prefix = to_wsl_path(args.ml_prefix)
        if args.support_prefix:
            args.support_prefix = to_wsl_path(args.support_prefix)
        log_file = to_wsl_path(log_file)
    # Step 1: ML tree search/model selection
    ml_prefix = args.ml_prefix or os.path.splitext(args.output)[0]
    ml_prefix = re.sub(r'_support$', '', ml_prefix)  # Remove _support if present
    ml_treefile = f"{ml_prefix}.treefile"
    ml_iqtree = f"{ml_prefix}.iqtree"
    ml_cmd = [
        "iqtree2",
        "-s", args.alignment,
        "-m", "MFP",
        "-nt", "AUTO",
        "-n", str(args.n),
        "--ninit", str(args.ninit),
        "--safe",
        "--seed", str(args.seed),
        "--keep-ident",
        "--redo",
        "-pre", ml_prefix
    ]
    with open(log_file, 'w') as log:
        log.write(f"Step 1: ML tree search/model selection\n")
    run_cmd(ml_cmd, log_file, use_wsl=use_wsl)
    # Step 2: Support estimation
    best_model = extract_best_fit_model(ml_iqtree)
    support_prefix = args.support_prefix or (ml_prefix + "_support")
    support_cmd = [
        "iqtree2",
        "-s", args.alignment,
        "-m", best_model,
        "-nt", "AUTO",
        "--keep-ident",
        "--safe",
        "--seed", str(args.seed),
        "-bb", str(args.bb),
        "-alrt", str(args.alrt),
        "--redo",
        "-pre", support_prefix
    ]
    # Always add --nstop (default 100 if not provided)
    nstop_val = args.nstop if args.nstop is not None else 100
    support_cmd.extend(["--nstop", str(nstop_val)])
    with open(log_file, 'a') as log:
        log.write(f"\nStep 2: Branch support estimation\n")
    run_cmd(support_cmd, log_file, use_wsl=use_wsl)
    # Copy final support tree to requested output
    final_tree = f"{support_prefix}.treefile"
    # Only copy if source and destination are different
    if os.path.abspath(final_tree) != os.path.abspath(args.output):
        import shutil
        shutil.copy(final_tree, args.output)
    else:
        with open(log_file, 'a') as log:
            log.write(f"Output tree file {final_tree} is already at the requested location. No copy needed.\n")
        open(args.output, 'w').close()

if __name__ == "__main__":
    main()
