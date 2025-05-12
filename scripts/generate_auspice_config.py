#!/usr/bin/env python3
"""
Dynamically generate or modify Auspice configuration JSON based on metadata and pathogen specifics.
This script enhances the static auspice_config.json with dynamic elements based on the dataset.
"""

import os
import sys
import json
import argparse
import logging
import pandas as pd
import yaml
from pathlib import Path
from collections import Counter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("generate_auspice_config")

def load_yaml(file_path):
    """Load YAML file into Python dictionary."""
    try:
        with open(file_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading YAML file {file_path}: {e}")
        return None

def load_json(file_path):
    """Load JSON file into Python dictionary."""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error loading JSON file {file_path}: {e}")
        return None

def load_metadata(file_path):
    """Load metadata TSV file into pandas DataFrame."""
    try:
        return pd.read_csv(file_path, sep='\t')
    except Exception as e:
        logger.error(f"Error loading metadata file {file_path}: {e}")
        return None

def get_unique_values(df, column, min_count=1):
    """Get unique values from a column with at least min_count occurrences."""
    if column not in df.columns:
        return []
    
    # Count occurrences of each value
    value_counts = Counter(df[column].dropna())
    
    # Filter for values with at least min_count occurrences
    return [val for val, count in value_counts.items() if count >= min_count]

def derive_color_scale(values, color_scheme="tableau10"):
    """Generate a color scale for categorical values."""
    import colorcet as cc
    
    # Simple predefined color schemes
    schemes = {
        "tableau10": [
            "#4e79a7", "#f28e2c", "#e15759", "#76b7b2", 
            "#59a14f", "#edc949", "#af7aa1", "#ff9da7", 
            "#9c755f", "#bab0ab"
        ],
        "category10": [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
            "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
            "#bcbd22", "#17becf"
        ]
    }
    
    # If scheme is a named colorcet palette, use that
    if color_scheme.startswith("cc_") and color_scheme[3:] in dir(cc):
        palette = getattr(cc, color_scheme[3:])
        # Convert to hex codes
        colors = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b, _ in palette]
    else:
        colors = schemes.get(color_scheme, schemes["tableau10"])
    
    # Cycle colors if more values than colors
    result = []
    for i, value in enumerate(values):
        color_idx = i % len(colors)
        result.append([str(value), colors[color_idx]])
    
    return result

def get_categorical_metadata_fields(df, exclude_columns=None, min_unique=2, max_unique=50):
    """
    Identify suitable categorical metadata fields for coloring.
    
    Args:
        df: Metadata DataFrame
        exclude_columns: List of columns to exclude
        min_unique: Minimum number of unique values
        max_unique: Maximum number of unique values for a field to be considered categorical
        
    Returns:
        List of column names that are good candidates for categorical coloring
    """
    if exclude_columns is None:
        exclude_columns = ["strain", "accession", "gisaid_epi_isl", "sequence_name"]
    
    categorical_fields = []
    
    for column in df.columns:
        if column in exclude_columns:
            continue
        
        # Skip columns with too many missing values
        missing_ratio = df[column].isna().mean()
        if missing_ratio > 0.7:  # Skip if more than 70% missing
            continue
        
        unique_values = df[column].nunique()
        if min_unique <= unique_values <= max_unique:
            categorical_fields.append(column)
    
    return categorical_fields

def enhance_auspice_config(base_config, metadata_df, pathogen_config):
    """
    Enhance the base Auspice config with dynamic elements based on metadata.
    
    Args:
        base_config: Base Auspice configuration dictionary
        metadata_df: Pandas DataFrame with metadata
        pathogen_config: Pathogen-specific configuration
        
    Returns:
        Enhanced Auspice configuration dictionary
    """
    # Create a copy to avoid modifying the original
    config = base_config.copy()
    
    # If metadata is empty, return the base config unchanged
    if metadata_df is None or metadata_df.empty:
        logger.warning("No metadata available, using base config without enhancement")
        return config
    
    # Add title and maintainers if not present
    if "title" not in config:
        config["title"] = f"{pathogen_config.get('name', 'Pathogen')} Genomic Surveillance"
    
    # Enhance colorings if not already comprehensive
    if "colorings" not in config:
        config["colorings"] = []
    
    existing_coloring_keys = [c.get("key") for c in config["colorings"]]
    
    # Check for additional categorical fields that could be useful for coloring
    categorical_fields = get_categorical_metadata_fields(metadata_df)
    
    for field in categorical_fields:
        if field not in existing_coloring_keys:
            unique_values = get_unique_values(metadata_df, field)
            # Only add fields with at least 2 values
            if len(unique_values) >= 2:
                logger.info(f"Adding new coloring for field: {field}")
                
                # Generate color scale for values
                color_scale = derive_color_scale(unique_values)
                
                config["colorings"].append({
                    "key": field,
                    "title": field.replace("_", " ").title(),
                    "type": "categorical",
                    "scale": color_scale
                })
    
    # Ensure all geo_resolutions in metadata are available for mapping
    geo_fields = [col for col in metadata_df.columns if col in ["region", "country", "division", "location"]]
    existing_geo = config.get("geo_resolutions", [])
    
    # Add any missing geo resolutions
    for field in geo_fields:
        if field not in existing_geo:
            logger.info(f"Adding geo resolution: {field}")
            config.setdefault("geo_resolutions", []).append(field)
    
    # Include special handling for RVF segments if applicable
    has_segments = pathogen_config.get("has_segments", False)
    if has_segments and "segment" in metadata_df.columns:
        segments = get_unique_values(metadata_df, "segment")
        if segments and len(segments) > 1:
            # Ensure segment coloring is present and customized
            segment_coloring = next((c for c in config["colorings"] if c.get("key") == "segment"), None)
            if not segment_coloring:
                # Create segment coloring with custom colors
                segment_colors = []
                for segment in segments:
                    if segment == "L":
                        segment_colors.append(["L", "#8F2727"])
                    elif segment == "M":
                        segment_colors.append(["M", "#3F548F"])
                    elif segment == "S":
                        segment_colors.append(["S", "#458F3F"])
                    else:
                        segment_colors.append([segment, "#999999"])
                
                config["colorings"].append({
                    "key": "segment",
                    "title": "Genome Segment",
                    "type": "categorical",
                    "scale": segment_colors
                })
                
                # Set display defaults to highlight segments if not set
                if "display_defaults" not in config:
                    config["display_defaults"] = {
                        "color_by": "segment"
                    }
                elif "color_by" not in config["display_defaults"]:
                    config["display_defaults"]["color_by"] = "segment"
    
    # Add filters based on metadata columns
    if "filters" not in config:
        config["filters"] = []
        
    # Standard filters to include
    standard_filters = ["region", "country", "division", "host", "author", "segment"]
    
    # Add categorical fields as filters if not already present
    filter_candidates = set(categorical_fields + standard_filters)
    existing_filters = set(config["filters"])
    
    for filter_field in filter_candidates:
        if filter_field in metadata_df.columns and filter_field not in existing_filters:
            config["filters"].append(filter_field)
    
    # Ensure standard visualization panels
    if "panels" not in config:
        config["panels"] = ["tree", "map", "entropy"]
    
    # Add metadata descriptions if not present
    if "metadata" not in config:
        config["metadata"] = []
    
    # Get existing metadata fields
    existing_metadata_fields = [m.get("name") for m in config.get("metadata", [])]
    
    # Add descriptions for key fields if not already present
    key_fields = {
        "segment": "Genome segment of the virus",
        "host": "Host species from which the sample was collected",
        "country": "Country where the sample was collected",
        "region": "Geographic region of sample collection",
        "author": "Authors of the sequence submission",
        "clade_membership": "Phylogenetic clade assignment"
    }
    
    for field, description in key_fields.items():
        if field in metadata_df.columns and field not in existing_metadata_fields:
            config["metadata"].append({
                "name": field,
                "type": "categorical", 
                "description": description
            })
    
    # Add description if missing
    if "description" not in config:
        config["description"] = f"{pathogen_config.get('name', 'Pathogen')} genomic surveillance dashboard. This visualization enables tracking of genetic lineages across time and geography."
    
    return config

def main():
    parser = argparse.ArgumentParser(description="Generate enhanced Auspice configuration JSON")
    parser.add_argument("--base-config", required=True, help="Base Auspice config JSON file")
    parser.add_argument("--metadata", required=True, help="Metadata TSV file")
    parser.add_argument("--master-config", required=True, help="Master config YAML file")
    parser.add_argument("--pathogen", help="Pathogen ID (defaults to active_pathogen in master config)")
    parser.add_argument("--output", required=True, help="Output enhanced config JSON file")
    args = parser.parse_args()
    
    # Load master config
    master_config = load_yaml(args.master_config)
    if not master_config:
        sys.exit(1)
    
    # Determine which pathogen to use
    pathogen_id = args.pathogen if args.pathogen else master_config.get("active_pathogen")
    if not pathogen_id or pathogen_id not in master_config.get("pathogens", {}):
        logger.error(f"Invalid pathogen ID: {pathogen_id}")
        sys.exit(1)
    
    pathogen_config = master_config["pathogens"][pathogen_id]
    
    # Load base Auspice config
    base_config = load_json(args.base_config)
    if not base_config:
        sys.exit(1)
    
    # Load metadata
    metadata_df = load_metadata(args.metadata)
    if metadata_df is None:
        logger.warning("Failed to load metadata, proceeding with base config only")
        metadata_df = pd.DataFrame()
    
    # Generate enhanced config
    enhanced_config = enhance_auspice_config(base_config, metadata_df, pathogen_config)
    
    # Write enhanced config to output file
    try:
        with open(args.output, 'w') as f:
            json.dump(enhanced_config, f, indent=2)
        logger.info(f"Enhanced Auspice config successfully written to {args.output}")
    except Exception as e:
        logger.error(f"Error writing enhanced config to {args.output}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()