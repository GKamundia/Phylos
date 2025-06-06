#!/usr/bin/env python3
"""
Create Auspice JSON file directly from filtered sequences for visualization
"""

import os
import sys
import json
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Create Auspice JSON for RVF sequences")
    parser.add_argument("--output-dir", default="auspice", help="Output directory for Auspice JSON")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Define paths
    filtered_sequences = "results/filtered/rvf_filtered.fasta"
    filtered_metadata = "results/filtered/rvf_metadata.tsv"
    auspice_config = "config/auspice_config.json"
    output_json = os.path.join(args.output_dir, "rvf.json")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Check if input files exist
    if not os.path.exists(filtered_sequences):
        print(f"Error: Filtered sequences file {filtered_sequences} not found")
        sys.exit(1)
    
    if not os.path.exists(filtered_metadata):
        print(f"Error: Filtered metadata file {filtered_metadata} not found")
        sys.exit(1)
    
    # Create a simple Auspice config if it doesn't exist
    if not os.path.exists(auspice_config):
        print(f"Creating default Auspice config file at {auspice_config}")
        os.makedirs(os.path.dirname(auspice_config), exist_ok=True)
        
        default_config = {
            "title": "Rift Valley Fever Virus",
            "maintainers": [{"name": "Nextstrain Team"}],
            "build_url": "https://github.com/example/rvf-nextstrain",
            "colorings": [
                {"key": "Segment", "title": "Segment", "type": "categorical"},
                {"key": "country", "title": "Country", "type": "categorical"},
                {"key": "host", "title": "Host", "type": "categorical"},
                {"key": "date", "title": "Date", "type": "categorical"},
                {"key": "num_date", "title": "Date", "type": "continuous"}
            ],
            "geo_resolutions": [
                {"key": "country", "title": "Country"}
            ],
            "filters": ["country", "host", "Segment", "Nuc_Completeness"],
            "panels": ["tree", "map", "entropy"]
        }
        
        with open(auspice_config, 'w') as f:
            json.dump(default_config, f, indent=2)
    
    # Create a basic tree JSON for Auspice
    print("Creating Auspice JSON file...")
    
    try:
        # Read metadata
        metadata = pd.read_csv(filtered_metadata, sep='\t')
        
        # Get sequences
        sequences = list(SeqIO.parse(filtered_sequences, "fasta"))
          # Create a simple tree structure
        auspice_json = {
            "version": "v2",
            "meta": {
                "title": "Rift Valley Fever Virus",
                "updated": pd.Timestamp.now().strftime("%Y-%m-%d"),
                "description": "Rift Valley Fever Virus sequences",
                "colorings": [
                    {"key": "Segment", "title": "Segment", "type": "categorical"},
                    {"key": "Country", "title": "Country", "type": "categorical"},
                    {"key": "Host", "title": "Host", "type": "categorical"},
                    {"key": "Collection_Date", "title": "Collection Date", "type": "categorical"},
                    {"key": "Nuc_Completeness", "title": "Completeness", "type": "categorical"}
                ],
                "filters": ["Country", "Host", "Segment", "Nuc_Completeness"],
                "geo_resolutions": [
                    {"key": "country", "title": "Country", "demes": {}}
                ]
            },
            "tree": {
                "name": "ROOT",
                "children": []
            },
            "nodes": {}
        }          
        # Add countries to geo_resolutions
        if 'Country' in metadata.columns:
            countries = metadata['Country'].dropna().unique()
            for country in countries:
                if country and str(country).strip():
                    auspice_json["meta"]["geo_resolutions"][0]["demes"][country] = {"latitude": 0, "longitude": 0}
        
        # Create a flat tree (all nodes connected to root)
        for i, seq in enumerate(sequences):
            seq_id = seq.id
            node_name = f"NODE_{i}"
            
            # Get metadata for this sequence
            seq_meta = metadata[metadata['strain'] == seq_id]
            
            if not seq_meta.empty:
                node_attrs = {}
                
                # Add segment info
                if 'Segment' in seq_meta.columns:
                    segment = seq_meta['Segment'].values[0]
                    node_attrs["Segment"] = {"value": segment}
                
                # Add country info
                if 'Country' in seq_meta.columns:
                    country = seq_meta['Country'].values[0]
                    if not pd.isna(country) and str(country).strip():
                        node_attrs["country"] = {"value": str(country)}
                
                # Add host info
                if 'Host' in seq_meta.columns:
                    host = seq_meta['Host'].values[0]
                    if not pd.isna(host) and str(host).strip():
                        node_attrs["host"] = {"value": str(host)}
                
                # Add date info
                if 'Collection_Date' in seq_meta.columns:
                    date = seq_meta['Collection_Date'].values[0]
                    if not pd.isna(date) and str(date).strip():
                        node_attrs["date"] = {"value": str(date)}
                
                # Add sequence info
                node_attrs["sequence"] = {"value": str(seq.seq)}
                
                # Add node to tree
                auspice_json["tree"]["children"].append({
                    "name": node_name,
                    "branch_length": 0.0001
                })
                
                # Add node attributes
                auspice_json["nodes"][node_name] = {
                    "name": seq_id,
                    "node_attrs": node_attrs
                }
        
        # Write Auspice JSON
        with open(output_json, 'w') as f:
            json.dump(auspice_json, f, indent=2)
        
        print(f"Auspice JSON file created at {output_json}")
        print(f"Use 'auspice view --datasetDir {args.output_dir}' to visualize")
        
    except Exception as e:
        print(f"Error creating Auspice JSON: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
