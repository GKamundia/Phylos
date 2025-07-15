#!/usr/bin/env python3
"""
Compare tip names in a Newick tree with strain names in a metadata TSV file.
Reports mismatches to help ensure all tree tips have corresponding metadata.

Usage:
  python compare_tree_and_metadata.py \
    --tree results/segments/M/refined/rvf_M_refined.tree \
    --metadata results/segments/M/filtered/rvf_M_metadata_latlong.tsv
"""
import argparse
from collections import Counter

try:
    from Bio import Phylo
except ImportError:
    print("Biopython is required. Install with: pip install biopython")
    exit(1)
import csv

def get_tree_tips(tree_path):
    tree = Phylo.read(tree_path, "newick")
    return set([term.name for term in tree.get_terminals()])

def get_metadata_strains(metadata_path):
    with open(metadata_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        strains = [row['strain'] for row in reader if 'strain' in row]
    return set(strains)

def main(tree_path, metadata_path):
    tree_tips = get_tree_tips(tree_path)
    metadata_strains = get_metadata_strains(metadata_path)

    tips_missing_in_metadata = tree_tips - metadata_strains
    strains_missing_in_tree = metadata_strains - tree_tips

    print(f"Total tips in tree: {len(tree_tips)}")
    print(f"Total strains in metadata: {len(metadata_strains)}\n")

    print(f"Tips in tree but missing from metadata ({len(tips_missing_in_metadata)}):")
    for tip in sorted(tips_missing_in_metadata):
        print(f"  {tip}")
    print()

    print(f"Strains in metadata but missing from tree ({len(strains_missing_in_tree)}):")
    for strain in sorted(strains_missing_in_tree):
        print(f"  {strain}")
    print()

    if not tips_missing_in_metadata and not strains_missing_in_tree:
        print("All tree tips have matching metadata strains. Ready for augur traits.")
    else:
        print("\nResolve these mismatches before re-running augur traits.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare tree tips and metadata strains.")
    parser.add_argument('--tree', required=True, help='Path to Newick tree file')
    parser.add_argument('--metadata', required=True, help='Path to metadata TSV file (with strain column)')
    args = parser.parse_args()
    main(args.tree, args.metadata)
