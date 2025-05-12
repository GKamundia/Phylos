#!/usr/bin/env python3
"""
Temporarily set the active pathogen in master_config.yaml
"""

import os
import sys
import shutil
import yaml

def set_active_pathogen(pathogen_name):
    """Set the active_pathogen value in master_config.yaml and create a backup."""
    config_path = "config/master_config.yaml"
    backup_path = "config/master_config.yaml.bak"
    
    # First check if pathogen exists in config
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        if pathogen_name not in config.get("pathogens", {}):
            print(f"Error: Pathogen '{pathogen_name}' not found in master_config.yaml")
            sys.exit(1)
        
        # Create a backup if it doesn't exist
        if not os.path.exists(backup_path):
            shutil.copy2(config_path, backup_path)
            print(f"Created backup of original config at {backup_path}")
        
        # Update active_pathogen
        config["active_pathogen"] = pathogen_name
        
        # Write updated config
        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        
        print(f"Set active_pathogen to '{pathogen_name}'")
        
    except Exception as e:
        print(f"Error modifying config: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python set_active_pathogen.py PATHOGEN_NAME")
        sys.exit(1)
    
    set_active_pathogen(sys.argv[1])