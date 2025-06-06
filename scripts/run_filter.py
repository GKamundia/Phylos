#!/usr/bin/env python3
"""
Run augur filter with error handling
"""

import subprocess
import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--metadata", required=True) 
    parser.add_argument("--output-sequences", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(os.path.dirname(args.output_sequences), exist_ok=True)
    
    # Run augur filter
    cmd = [
        "augur", "filter",
        "--sequences", args.sequences,
        "--metadata", args.metadata,
        "--output-sequences", args.output_sequences,
        "--min-length", "4",
        "--exclude-where", "Host=''",
        "--exclude-where", "Collection_Date=''", 
        "--include-where", "Nuc_Completeness='complete'",
        "--output-metadata", args.output_metadata
    ]
    
    with open(args.log, 'w') as log_file:
        log_file.write("Filtering for complete sequences across all segments\n")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            log_file.write(result.stdout)
            log_file.write(result.stderr)
            
            if result.returncode != 0:
                print(f"Augur filter failed with return code {result.returncode}")
                return result.returncode
                
        except Exception as e:
            log_file.write(f"Error running augur filter: {e}\n")
            return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
