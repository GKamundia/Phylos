"""
Master Snakefile for Pathogen-Agnostic Nextstrain Pipeline
"""

import os
import sys
import yaml
from pathlib import Path
from snakemake.utils import min_version

# Ensure minimum Snakemake version for compatibility
min_version("7.18.2")

# ==== Configuration Loading ====

def load_yaml(yaml_file):
    """Load YAML file with error handling"""
    try:
        with open(yaml_file) as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading {yaml_file}: {e}", file=sys.stderr)
        sys.exit(1)

# Load and merge configurations
master_config = load_yaml("config/master_config.yaml")
active_pathogen = master_config["active_pathogen"]

try:
    pathogen_config = master_config["pathogens"][active_pathogen]
except KeyError:
    print(f"Error: Pathogen '{active_pathogen}' not found in master_config.yaml", file=sys.stderr)
    sys.exit(1)

# Load pathogen-specific configuration
pathogen_config_path = pathogen_config["config_path"]
if not os.path.exists(pathogen_config_path):
    print(f"Error: Pathogen config file {pathogen_config_path} not found", file=sys.stderr)
    sys.exit(1)

pathogen_specific_config = load_yaml(pathogen_config_path)

# Merge configurations with proper fallbacks and inheritance
config = {
    "pathogen": active_pathogen,
    "pathogen_name": pathogen_config["name"],
    "has_segments": pathogen_config["has_segments"],
    "data": pathogen_specific_config.get("data", {}),
    "filter": {**master_config["common"]["filter"], **pathogen_specific_config.get("filter", {})},
    "subsample": {**master_config["common"]["subsample"], **pathogen_specific_config.get("subsample", {})},
    "refine": {**master_config["common"].get("refine", {}), **pathogen_specific_config.get("refine", {})},
    "update": master_config["common"]["update"],
    "qc": master_config["common"].get("qc", {}),
    "resources": {**master_config["common"].get("resources", {}), **pathogen_specific_config.get("resources", {})},
    "workflow": {**master_config["common"].get("workflow", {}), **pathogen_specific_config.get("workflow", {})},
    "export": master_config["common"].get("export", {}),  # Add export configuration
}

# ==== Workflow Configuration ====

# Define output prefix based on active pathogen
output_prefix = active_pathogen

# Handle segment configuration with enhanced conditionality
segments = []
if config["has_segments"]:
    segment_mode = config["data"].get("segment_mode", "single")
    specified_segment = config["data"].get("segment", "")
    
    # Determine segments based on configuration
    if segment_mode == "single" and specified_segment:
        segments = [specified_segment]
    elif segment_mode == "multi" or specified_segment == "all":
        # Get segments from pathogen-specific config or use default
        segments = pathogen_specific_config.get("segments", ["L", "M", "S"])
        segment_mode = "multi"
    else:
        # Default to L segment if nothing specified
        segments = ["L"]
        segment_mode = "single"
else:
    # Non-segmented pathogen
    segments = [""]
    segment_mode = "single"

# Export for use in rule files
config["segment_mode"] = segment_mode
config["segments"] = segments

# ==== Directory Setup ====

# Create required directories
required_dirs = [
    "logs", 
    "results",
    Path("data/sequences/raw"),
    Path("data/metadata/raw"),
    "results/nextclade",
    "results/qc_reports",
    "benchmarks"
]

# Create segment-specific directories if needed
if segment_mode == "multi":
    for segment in segments:
        required_dirs.extend([
            f"results/segments/{segment}",
            f"results/segments/{segment}/raw",
            f"results/segments/{segment}/filtered", 
            f"results/segments/{segment}/aligned",
            f"results/segments/{segment}/tree",
            f"results/segments/{segment}/node_data",
            f"results/segments/{segment}/nextclade",
            f"results/segments/{segment}/qc_reports",
            f"results/segments/{segment}/auspice"
        ])

# Create all directories
for directory in required_dirs:
    os.makedirs(directory, exist_ok=True)

# ==== Include Rule Modules ====

# Include modular rule files based on configuration
include: "workflow/core_rules/data_acquisition.smk"
include: "workflow/core_rules/filtering.smk"
include: "workflow/core_rules/alignment.smk"
include: "workflow/core_rules/analysis.smk"
include: "workflow/core_rules/export.smk"
include: "workflow/utils/checkpoint_handlers.smk"
include: "workflow/utils/validation.smk"
include: "workflow/utils/notification.smk"

# Include segment-specific rule modules conditionally
if segment_mode == "multi":
    include: "workflow/segment_rules/multi_segment.smk"
else:
    include: "workflow/segment_rules/single_segment.smk"

# ==== Output Target Rules ====

# Define get_final_outputs function BEFORE rule all
def get_final_outputs(wildcards):
    """Generate a list of final outputs based on configuration"""
    outputs = []
    
    # QC reports always included
    outputs.append(f"results/qc_reports/{output_prefix}_qc_summary.json")
    
    # Single segment mode outputs
    if segment_mode == "single":
        outputs.append(f"results/auspice/{output_prefix}.json")
    
    # Multi-segment mode outputs
    else:
        # Individual segment outputs
        for segment in segments:
            outputs.append(f"results/segments/{segment}/auspice/{output_prefix}_{segment}.json")
        
        # Combined visualization if enabled
        if config.get("workflow", {}).get("create_combined_view", True):
            outputs.append(f"results/auspice/{output_prefix}_combined.json")
    
    return outputs

rule all:
    input:
        # Dynamic output targets based on segment mode
        get_final_outputs

# Define a target rule for validated outputs
rule validate_outputs:
    input:
        get_final_outputs
    output:
        touch("results/.validation_complete")
    run:
        print(f"Pipeline completed successfully for {config['pathogen_name']} ({config['pathogen']})")
        print(f"Segment mode: {segment_mode}")
        if segment_mode == "multi":
            print(f"Analyzed segments: {', '.join(segments)}")
        print("\nResults available at:")
        for file in input:
            print(f"  - {file}")