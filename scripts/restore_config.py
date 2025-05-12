#!/usr/bin/env python3
"""
Restore master_config.yaml from backup after a scheduled run
"""

import os
import sys
import shutil

def restore_config():
    """Restore master_config.yaml from backup if it exists."""
    config_path = "config/master_config.yaml"
    backup_path = "config/master_config.yaml.bak"
    
    if not os.path.exists(backup_path):
        print("No backup file found, nothing to restore")
        return
    
    try:
        shutil.copy2(backup_path, config_path)
        os.remove(backup_path)
        print(f"Restored {config_path} from backup and removed {backup_path}")
    except Exception as e:
        print(f"Error restoring config: {e}")
        sys.exit(1)

if __name__ == "__main__":
    restore_config()