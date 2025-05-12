#!/usr/bin/env python3
"""
Combine Auspice JSONs from multiple segments into a single visualization
"""

import json
import logging

import snakemake

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get input parameters from Snakemake
segment_jsons = snakemake.input.segment_jsons
combined_json = snakemake.output.combined_json
title = snakemake.params.title
log_file = snakemake.log[0]

# Configure logging to file
file_handler = logging.FileHandler(log_file)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(file_handler)

# Function to load JSON file
def load_json(filename):
    with open(filename, 'r') as f:
        return json.load(f)

# Load segment JSONs
logger.info(f"Loading segment JSONs: {segment_jsons}")
segment_data = [load_json(f) for f in segment_jsons]

# Extract segment names from filenames
import os
segment_names = [os.path.basename(f).split('.')[0].split('_')[-1] for f in segment_jsons]
logger.info(f"Detected segments: {segment_names}")

# Start with the first segment's data as the template
combined = segment_data[0].copy()

# Update metadata
combined["meta"]["title"] = title
combined["meta"]["segments"] = segment_names

# Create a dataset section for each segment
combined["datasets"] = {}
for i, (name, data) in enumerate(zip(segment_names, segment_data)):
    # Extract key elements from each segment
    combined["datasets"][name] = {
        "name": f"Segment {name}",
        "description": data["meta"].get("description", f"Analysis of {name} segment"),
        "tree": data["tree"],
        "entropy": data.get("entropy", {}),
        "frequencies": data.get("frequencies", {}),
    }

# Write combined JSON
logger.info(f"Writing combined JSON to {combined_json}")
with open(combined_json, 'w') as f:
    json.dump(combined, f, indent=2)

logger.info("Combined JSON created successfully")