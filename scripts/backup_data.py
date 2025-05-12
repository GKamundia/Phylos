#!/usr/bin/env python3
"""
Automated backup script for critical RVF-Nextstrain data and configurations
"""

import os
import sys
import shutil
import logging
import argparse
import datetime
import json
import yaml
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("backup_data")

# Default backup configuration
DEFAULT_CONFIG = {
    "critical_data": [
        "data/sequences/raw",
        "data/metadata/raw",
        "data/metadata/*.tsv"
    ],
    "critical_configs": [
        "config/*.yaml",
        "config/*.json",
        "config/*.tsv",
        "pathogens/*/config/*.yaml",
        "pathogens/*/config/*.json",
        "pathogens/*/config/*.tsv"
    ],
    "output_dirs": [
        "results/auspice",
        "results/qc_reports"
    ],
    "backup_dir": "backups",
    "retention_days": 30,
    "retention_count": 10,
    "compression": True,
    "remote_sync": False,
    "remote_path": None
}

def read_master_config():
    """Read the master configuration file to get pathogen info"""
    try:
        config_path = "config/master_config.yaml"
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.warning(f"Could not read master config: {e}")
        return {}

def create_backup(config_file=None, backup_name=None):
    """Create a backup of critical data and configurations"""
    # Load backup configuration
    if config_file and os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                config = json.load(f)
        except Exception as e:
            logger.error(f"Error loading backup config: {e}")
            config = DEFAULT_CONFIG
    else:
        config = DEFAULT_CONFIG

    # Create backup directory with timestamp
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = Path(config["backup_dir"])
    
    if backup_name:
        backup_path = backup_dir / f"{backup_name}_{timestamp}"
    else:
        # Read master config to get active pathogen
        master_config = read_master_config()
        active_pathogen = master_config.get("active_pathogen", "unknown")
        backup_path = backup_dir / f"{active_pathogen}_backup_{timestamp}"
    
    backup_path.mkdir(parents=True, exist_ok=True)
    logger.info(f"Creating backup in: {backup_path}")
    
    # Create metadata about the backup
    backup_metadata = {
        "timestamp": timestamp,
        "created": datetime.datetime.now().isoformat(),
        "pathogen": master_config.get("active_pathogen", "unknown"),
        "content": []
    }
    
    # Backup critical data
    for data_pattern in config["critical_data"]:
        # Use Path.glob for pattern matching
        for path in Path().glob(data_pattern):
            if path.is_file():
                # For files, create directory structure and copy file
                target_path = backup_path / path
                target_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target_path)
                logger.info(f"Backed up file: {path}")
                backup_metadata["content"].append({"type": "file", "path": str(path)})
            elif path.is_dir():
                # For directories, copy the entire directory tree
                target_path = backup_path / path
                target_path.parent.mkdir(parents=True, exist_ok=True)
                if os.path.exists(path):
                    try:
                        shutil.copytree(path, target_path)
                        logger.info(f"Backed up directory: {path}")
                        backup_metadata["content"].append({"type": "directory", "path": str(path)})
                    except Exception as e:
                        logger.error(f"Failed to backup {path}: {e}")
    
    # Backup configuration files
    config_backup_dir = backup_path / "config"
    config_backup_dir.mkdir(parents=True, exist_ok=True)
    
    for config_pattern in config["critical_configs"]:
        for config_file in Path().glob(config_pattern):
            if config_file.is_file():
                # Preserve directory structure under config/
                rel_path = config_file.relative_to(config_file.parts[0])
                target_file = config_backup_dir / rel_path
                target_file.parent.mkdir(parents=True, exist_ok=True)
                try:
                    shutil.copy2(config_file, target_file)
                    logger.info(f"Backed up config: {config_file}")
                    backup_metadata["content"].append({"type": "config", "path": str(config_file)})
                except Exception as e:
                    logger.error(f"Failed to backup {config_file}: {e}")
    
    # Backup selected output directories if requested
    if config.get("backup_outputs", False):
        for output_pattern in config.get("output_dirs", []):
            for output_path in Path().glob(output_pattern):
                if output_path.is_dir():
                    target_path = backup_path / output_path
                    target_path.parent.mkdir(parents=True, exist_ok=True)
                    try:
                        shutil.copytree(output_path, target_path)
                        logger.info(f"Backed up output: {output_path}")
                        backup_metadata["content"].append({"type": "output", "path": str(output_path)})
                    except Exception as e:
                        logger.error(f"Failed to backup output {output_path}: {e}")
    
    # Save backup metadata
    with open(backup_path / "backup_metadata.json", 'w') as f:
        json.dump(backup_metadata, f, indent=2)
    
    # Compress backup if configured
    if config.get("compression", False):
        try:
            archive_path = f"{backup_path}.tar.gz"
            shutil.make_archive(
                str(backup_path), 
                'gztar', 
                root_dir=str(backup_dir),
                base_dir=backup_path.name
            )
            logger.info(f"Created compressed archive: {archive_path}")
            
            # Remove uncompressed directory if compression was successful
            shutil.rmtree(backup_path)
            backup_path = Path(f"{archive_path}")
        except Exception as e:
            logger.error(f"Failed to compress backup: {e}")
    
    # Sync to remote location if configured
    if config.get("remote_sync", False) and config.get("remote_path"):
        try:
            # Use rsync or similar to send to remote location
            remote_path = config["remote_path"]
            os.system(f"rsync -av {backup_path} {remote_path}")
            logger.info(f"Synced backup to remote location: {remote_path}")
        except Exception as e:
            logger.error(f"Failed to sync to remote location: {e}")
    
    # Clean up old backups based on retention policy
    clean_old_backups(config)
    
    logger.info(f"Backup completed: {backup_path}")
    return str(backup_path)

def clean_old_backups(config):
    """Clean up old backups based on retention policy"""
    backup_dir = Path(config["backup_dir"])
    if not backup_dir.exists():
        return
    
    # Get all backups ordered by creation date (oldest first)
    backups = []
    for item in backup_dir.glob("*_backup_*"):
        if item.is_dir() or (item.is_file() and item.suffix in ['.gz', '.zip']):
            try:
                backups.append((item, item.stat().st_mtime))
            except Exception:
                continue
    
    # Sort by modification time (oldest first)
    backups.sort(key=lambda x: x[1])
    
    # Remove old backups based on retention count
    retention_count = config.get("retention_count", 10)
    if len(backups) > retention_count:
        for backup, _ in backups[:-retention_count]:
            try:
                if backup.is_dir():
                    shutil.rmtree(backup)
                else:
                    backup.unlink()
                logger.info(f"Removed old backup: {backup}")
            except Exception as e:
                logger.error(f"Failed to remove old backup {backup}: {e}")
    
    # Remove backups older than retention_days
    retention_days = config.get("retention_days", 30)
    if retention_days > 0:
        cutoff_time = datetime.datetime.now() - datetime.timedelta(days=retention_days)
        cutoff_timestamp = cutoff_time.timestamp()
        
        for backup, mtime in backups:
            if mtime < cutoff_timestamp:
                try:
                    if backup.is_dir():
                        shutil.rmtree(backup)
                    else:
                        backup.unlink()
                    logger.info(f"Removed expired backup: {backup}")
                except Exception as e:
                    logger.error(f"Failed to remove expired backup {backup}: {e}")

def restore_backup(backup_path, target_dir=None, selective=False, data_only=False, config_only=False):
    """Restore data from a backup"""
    backup_path = Path(backup_path)
    
    # Handle compressed backups
    is_compressed = False
    if not backup_path.is_dir() and backup_path.suffix in ['.gz', '.zip']:
        is_compressed = True
        # Create a temporary extraction directory
        extract_dir = Path("backups/tmp_extract")
        extract_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            if backup_path.suffix == '.gz':
                shutil.unpack_archive(backup_path, extract_dir)
                # Find the extracted directory
                extracted_dirs = list(extract_dir.glob("*_backup_*"))
                if not extracted_dirs:
                    raise ValueError("Could not find backup directory in archive")
                backup_path = extracted_dirs[0]
            else:
                raise ValueError(f"Unsupported archive format: {backup_path.suffix}")
        except Exception as e:
            logger.error(f"Failed to extract backup: {e}")
            return False
    
    # Load backup metadata
    metadata_file = backup_path / "backup_metadata.json"
    try:
        if metadata_file.exists():
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
        else:
            metadata = {"content": []}
            # Try to infer content from directory structure
            for path in backup_path.rglob("*"):
                if path.is_file() and path.name != "backup_metadata.json":
                    rel_path = path.relative_to(backup_path)
                    # Skip metadata files
                    if str(rel_path) == "backup_metadata.json":
                        continue
                    parent_dir = str(rel_path.parts[0]) if rel_path.parts else ""
                    if parent_dir == "config":
                        metadata["content"].append({"type": "config", "path": str(rel_path)})
                    else:
                        metadata["content"].append({"type": "file", "path": str(rel_path)})
    except Exception as e:
        logger.error(f"Failed to load backup metadata: {e}")
        metadata = {"content": []}
    
    # Handle selective restore
    if selective:
        # Here you'd implement an interactive selection process
        # For this example, we'll just proceed with all files
        logger.info("Selective restore requested, but not implemented. Restoring all files.")
    
    # Restore files based on filtering options
    for item in metadata.get("content", []):
        item_type = item.get("type", "")
        item_path = item.get("path", "")
        
        # Skip based on filter options
        if data_only and item_type == "config":
            continue
        if config_only and item_type != "config":
            continue
            
        # Source path in backup
        if item_type == "config":
            source_path = backup_path / "config" / Path(item_path).relative_to(Path(item_path).parts[0])
        else:
            source_path = backup_path / item_path
        
        # Target path for restoration
        if target_dir:
            if item_type == "config":
                # For configs, maintain their structure
                target_path = Path(target_dir) / item_path
            else:
                # For data, restore into the target directory
                target_path = Path(target_dir) / Path(item_path).name
        else:
            # Restore to original location
            target_path = Path(item_path)
        
        # Create parent directories
        target_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Perform the restore
        try:
            if source_path.is_file():
                shutil.copy2(source_path, target_path)
                logger.info(f"Restored: {item_path} to {target_path}")
            elif source_path.is_dir():
                if target_path.exists():
                    shutil.rmtree(target_path)
                shutil.copytree(source_path, target_path)
                logger.info(f"Restored directory: {item_path} to {target_path}")
        except Exception as e:
            logger.error(f"Failed to restore {item_path}: {e}")
    
    # Clean up temporary extraction directory if necessary
    if is_compressed and extract_dir.exists():
        shutil.rmtree(extract_dir)
    
    logger.info("Restore completed")
    return True

def list_backups(verbose=False):
    """List available backups"""
    backup_dir = Path("backups")
    if not backup_dir.exists():
        logger.info("No backups found (backup directory does not exist)")
        return []
    
    backups = []
    for item in backup_dir.glob("*_backup_*"):
        if item.is_dir() or (item.is_file() and item.suffix in ['.gz', '.zip']):
            try:
                backup_info = {
                    "path": str(item),
                    "date": datetime.datetime.fromtimestamp(item.stat().st_mtime).isoformat(),
                    "size": item.stat().st_size,
                    "type": "directory" if item.is_dir() else "archive"
                }
                
                # Try to extract pathogen information from name
                name_parts = item.name.split("_backup_")
                if name_parts:
                    backup_info["pathogen"] = name_parts[0]
                
                # Load metadata for more details if verbose
                if verbose:
                    metadata_path = item / "backup_metadata.json" if item.is_dir() else None
                    if metadata_path and metadata_path.exists():
                        try:
                            with open(metadata_path, 'r') as f:
                                metadata = json.load(f)
                                backup_info["metadata"] = metadata
                        except Exception:
                            pass
                
                backups.append(backup_info)
                logger.info(f"Found backup: {item} ({backup_info['date']})")
            except Exception as e:
                logger.error(f"Error processing backup {item}: {e}")
    
    # Sort by date (newest first)
    backups.sort(key=lambda x: x["date"], reverse=True)
    return backups

def main():
    """Main function when running as a script"""
    parser = argparse.ArgumentParser(description="Backup and restore critical RVF-Nextstrain data")
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Backup command
    backup_parser = subparsers.add_parser("backup", help="Create a backup")
    backup_parser.add_argument("--config", help="Path to backup configuration file")
    backup_parser.add_argument("--name", help="Custom name for the backup")
    
    # Restore command
    restore_parser = subparsers.add_parser("restore", help="Restore from a backup")
    restore_parser.add_argument("backup", help="Backup path to restore from")
    restore_parser.add_argument("--target", help="Target directory for restoration")
    restore_parser.add_argument("--selective", action="store_true", help="Selectively restore files")
    restore_parser.add_argument("--data-only", action="store_true", help="Restore only data files")
    restore_parser.add_argument("--config-only", action="store_true", help="Restore only configuration files")
    
    # List command
    list_parser = subparsers.add_parser("list", help="List available backups")
    list_parser.add_argument("--verbose", "-v", action="store_true", help="Show detailed information")
    
    args = parser.parse_args()
    
    if args.command == "backup":
        backup_path = create_backup(args.config, args.name)
        print(f"Backup created: {backup_path}")
    
    elif args.command == "restore":
        if not args.backup:
            print("Error: Backup path is required for restore")
            return 1
        
        success = restore_backup(
            args.backup,
            args.target,
            args.selective,
            args.data_only,
            args.config_only
        )
        
        if success:
            print("Restore completed successfully")
            return 0
        else:
            print("Restore failed")
            return 1
    
    elif args.command == "list":
        backups = list_backups(args.verbose)
        if not backups:
            print("No backups found")
        else:
            print(f"Found {len(backups)} backup(s):")
            for i, backup in enumerate(backups):
                print(f"  [{i+1}] {backup['pathogen']} - {backup['date']} ({backup['path']})")
    
    else:
        parser.print_help()

if __name__ == "__main__":
    sys.exit(main())