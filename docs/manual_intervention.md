# Manual Intervention Guide for RVF-Nextstrain Pipeline

This guide explains how to manually intervene in the workflow when specific steps fail or need to be regenerated.

## Understanding Snakemake's Target-Based Execution

The Nextstrain pipeline uses Snakemake for workflow management, which offers several features for manual intervention:

1. **Target-specific execution**: Run only specific parts of the pipeline
2. **Force re-execution**: Re-run specific steps even if outputs exist
3. **Continue after failure**: Resume from where a failed run stopped

## Common Intervention Scenarios

### Scenario 1: Data Download Failed

If the data download step fails due to network issues or NCBI downtime:

```bash
# Manually retry just the download step
nextstrain build . --configfile config/master_config.yaml data/sequences/raw/rvf_sequences.fasta data/metadata/raw/rvf_metadata.tsv

# Alternative: force rerun the download rule
nextstrain build . --configfile config/master_config.yaml --forcerun download_data
```

### Scenario 2: Metadata Preparation Failed

If metadata validation or preparation failed:

```bash
# Regenerate only the metadata files
nextstrain build . --configfile config/master_config.yaml data/metadata/rvf_metadata.tsv

# Force rerun even if the file exists
nextstrain build . --configfile config/master_config.yaml --forcerun prepare_metadata
```

### Scenario 3: Alignment or Tree Building Failed

For failures in computationally intensive steps:

```bash
# Regenerate only the alignment
nextstrain build . --configfile config/master_config.yaml results/aligned/rvf_aligned.fasta

# Force rebuild the phylogenetic tree
nextstrain build . --configfile config/master_config.yaml --forcerun tree
```

### Scenario 4: Segment-Specific Issues (Multi-Segment Mode)

When problems occur with specific genome segments:

```bash
# Regenerate results for just the L segment
nextstrain build . --configfile config/master_config.yaml results/segments/L/auspice/rvf_L.json

# Force rerun the M segment alignment
nextstrain build . --configfile config/master_config.yaml --forcerun "align_segment_M"
```

### Scenario 5: Final Auspice Export Failed

If the final visualization export step fails:

```bash
# Regenerate only the auspice JSON files
nextstrain build . --configfile config/master_config.yaml results/auspice/rvf.json

# Force rerun the export rule
nextstrain build . --configfile config/master_config.yaml --forcerun export
```

## Resuming a Failed Run

When a workflow fails, you can resume from the last successful step:

```bash
# Resume the run (Snakemake will determine what needs to be rerun)
nextstrain build . --configfile config/master_config.yaml --keep-going
```

## Using Checkpoints to Resume Complex Workflows

For workflows with checkpoints (like multi-segment analysis), use:

```bash
# Resume after a checkpoint-dependent failure
nextstrain build . --configfile config/master_config.yaml --rerun-incomplete
```

## Cleaning Intermediate Files

If you need to start fresh from a specific point:

```bash
# Remove specific intermediate files
rm results/filtered/rvf_filtered.fasta
rm results/aligned/rvf_aligned.fasta

# Then run to regenerate just those files and their dependencies
nextstrain build . --configfile config/master_config.yaml results/aligned/rvf_aligned.fasta
```

## Advanced: Debugging with Dry Runs

To see what would be executed without actually running the commands:

```bash
# Dry run to see what would be executed
nextstrain build . --configfile config/master_config.yaml --dryrun results/auspice/rvf.json
```

## Recovering from Data Corruption

If you encounter data corruption in intermediate files:

1. Remove the corrupt files: rm [corrupt file path]
2. Force Snakemake to regenerate those files: nextstrain build . --configfile [master_config.yaml](http://_vscodecontentref_/1) [output file that depends on the corrupt file]
3. Data Recovery Backup Script

To ensure data safety, let's add a simple backup utility:

```python
#!/usr/bin/env python3
"""
Utilities for backing up critical data and configurations
"""

import os
import shutil
import logging
import datetime
import argparse
from pathlib import Path

logger = logging.getLogger(__name__)

def backup_critical_data(
    data_dirs=None,
    config_files=None,
    backup_dir=None,
    timestamp_format="%Y%m%d_%H%M%S"
):
    """
    Backup critical data and configuration files

    Args:
        data_dirs: List of data directories to backup
        config_files: List of config files to backup
        backup_dir: Directory to store backups
        timestamp_format: Format for timestamp in backup directory name

    Returns:
        Path to backup directory
    """
    # Default paths if None provided
    if data_dirs is None:
        data_dirs = [
            "data/sequences/raw",
            "data/metadata/raw"
        ]

    if config_files is None:
        config_files = [
            "config/master_config.yaml",
            "config/auspice_config.json",
            "config/metadata_schema.json",
            "config/data_tracking.json",
            "config/lat_longs.tsv"
        ]

    if backup_dir is None:
        backup_dir = "backups"

    # Create backup directory with timestamp
    timestamp = datetime.datetime.now().strftime(timestamp_format)
    backup_path = Path(backup_dir) / f"backup_{timestamp}"
    backup_path.mkdir(parents=True, exist_ok=True)

    # Backup data directories
    for data_dir in data_dirs:
        if os.path.exists(data_dir):
            target_dir = backup_path / data_dir
            target_dir.parent.mkdir(parents=True, exist_ok=True)
            try:
                shutil.copytree(data_dir, target_dir)
                logger.info(f"Backed up {data_dir} to {target_dir}")
            except Exception as e:
                logger.error(f"Failed to backup {data_dir}: {str(e)}")

    # Backup config files
    config_backup_dir = backup_path / "config"
    config_backup_dir.mkdir(parents=True, exist_ok=True)

    for config_file in config_files:
        if os.path.exists(config_file):
            target_file = config_backup_dir / Path(config_file).name
            try:
                shutil.copy2(config_file, target_file)
                logger.info(f"Backed up {config_file} to {target_file}")
            except Exception as e:
                logger.error(f"Failed to backup {config_file}: {str(e)}")

    return backup_path

def main():
    parser = argparse.ArgumentParser(description="Backup critical data and configurations")
    parser.add_argument("--backup-dir", default="backups", help="Directory to store backups")
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    backup_path = backup_critical_data(backup_dir=args.backup_dir)
    print(f"Backup completed: {backup_path}")

if __name__ == "__main__":
    main()
```

4. Integration with Workflow Monitoring
   Update the monitoring to check for failed steps and suggest manual recovery:

```python
# Add the following function
def generate_recovery_suggestions(error_files):
    """Generate recovery suggestions based on error patterns in log files"""
    suggestions = []

    for log_file in error_files:
        # Extract rule name from log file
        rule_match = re.search(r"logs/([^/]+)_", log_file)
        if not rule_match:
            continue

        rule_name = rule_match.group(1)

        if rule_name == "download_data":
            suggestions.append(
                "Data download failed. Try running:\n"
                "nextstrain build . --configfile config/master_config.yaml --forcerun download_data"
            )
        elif rule_name == "prepare_metadata":
            suggestions.append(
                "Metadata preparation failed. Try running:\n"
                "nextstrain build . --configfile config/master_config.yaml --forcerun prepare_metadata"
            )
        elif rule_name == "align":
            suggestions.append(
                "Sequence alignment failed. Try running:\n"
                "nextstrain build . --configfile config/master_config.yaml --forcerun align"
            )
        elif rule_name == "tree":
            suggestions.append(
                "Phylogenetic tree construction failed. Try running:\n"
                "nextstrain build . --configfile config/master_config.yaml --forcerun tree"
            )
        elif "export" in rule_name:
            suggestions.append(
                "Auspice export failed. Try running:\n"
                "nextstrain build . --configfile config/master_config.yaml --forcerun export"
            )
        elif "_segment_" in rule_name:
            # Extract segment name
            segment_match = re.search(r"_segment_([^_]+)", rule_name)
            if segment_match:
                segment = segment_match.group(1)
                suggestions.append(
                    f"Processing for segment {segment} failed. Try running:\n"
                    f"nextstrain build . --configfile config/master_config.yaml --forcerun {rule_name}"
                )

    # Add general recovery suggestions
    if error_files:
        suggestions.append(
            "To resume the entire workflow from where it failed:\n"
            "nextstrain build . --configfile config/master_config.yaml --keep-going"
        )

    return suggestions

# Update the generate_status_summary function to include recovery suggestions
def generate_status_summary(results: Dict[str, Any]) -> Dict[str, Any]:
    # Add this line near the end of the function:
    if error_files:
        status["recovery_suggestions"] = generate_recovery_suggestions(error_files)
```
