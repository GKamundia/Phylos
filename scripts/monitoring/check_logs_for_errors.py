# filepath: c:\Users\Anarchy\Documents\Data_Science\NextStrain\rvf-nextstrain\scripts\monitoring\check_logs_for_errors.py
#!/usr/bin/env python3
"""
Check logs for errors and warnings
"""

import os
import sys
import json
import re
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any

import snakemake

# Add parent directory to path for importing utils
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.log_utils import setup_logger, log_with_context

# Configure logger
logger = setup_logger(
    name="check_logs",
    log_file=snakemake.log[0] if "snakemake" in globals() else "check_logs.log"
)

def check_logs_for_patterns(
    log_files: List[str],
    error_patterns: List[str],
    warning_patterns: List[str]
) -> Dict[str, Any]:
    """
    Check log files for specified patterns
    
    Args:
        log_files: List of log file paths to check
        error_patterns: List of regex patterns that indicate errors
        warning_patterns: List of regex patterns that indicate warnings
        
    Returns:
        Dictionary with results
    """
    results = {
        "errors": [],
        "warnings": [],
        "status": "ok",
        "total_errors": 0,
        "total_warnings": 0,
        "checked_files": len(log_files),
        "timestamp": datetime.now().isoformat()
    }
    
    # Combine patterns for single-pass checking
    error_regex = re.compile("|".join(error_patterns), re.IGNORECASE)
    warning_regex = re.compile("|".join(warning_patterns), re.IGNORECASE)
    
    # Check each log file
    for log_file in log_files:
        try:
            file_errors = []
            file_warnings = []
            
            # Skip empty files
            if os.path.getsize(log_file) == 0:
                continue
                
            # Read log file
            with open(log_file, 'r', errors='replace') as f:
                line_number = 0
                for line in f:
                    line_number += 1
                    
                    # Check for errors
                    if error_regex.search(line):
                        # Get context (truncate if too long)
                        context = line.strip()
                        if len(context) > 200:
                            context = context[:197] + "..."
                            
                        file_errors.append({
                            "line": line_number,
                            "context": context
                        })
                    
                    # Check for warnings
                    elif warning_regex.search(line):
                        # Get context (truncate if too long)
                        context = line.strip()
                        if len(context) > 200:
                            context = context[:197] + "..."
                            
                        file_warnings.append({
                            "line": line_number,
                            "context": context
                        })
            
            # If errors or warnings found, add to results
            if file_errors:
                results["errors"].append({
                    "file": log_file,
                    "error_count": len(file_errors),
                    "errors": file_errors[:10]  # Limit to first 10 errors
                })
                results["total_errors"] += len(file_errors)
            
            if file_warnings:
                results["warnings"].append({
                    "file": log_file,
                    "warning_count": len(file_warnings),
                    "warnings": file_warnings[:10]  # Limit to first 10 warnings
                })
                results["total_warnings"] += len(file_warnings)
        
        except Exception as e:
            logger.error(f"Error checking log file {log_file}: {e}")
    
    # Set overall status
    if results["total_errors"] > 0:
        results["status"] = "error"
    elif results["total_warnings"] > 0:
        results["status"] = "warning"
    
    return results

def generate_status_summary(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate a simplified status summary
    
    Args:
        results: Results from check_logs_for_patterns
        
    Returns:
        Status summary dictionary
    """
    status = {
        "timestamp": results["timestamp"],
        "status": results["status"],
        "error_count": results["total_errors"],
        "warning_count": results["total_warnings"],
        "checked_files": results["checked_files"]
    }
    
    # Add details about most recent errors/warnings
    if results["errors"]:
        error_files = [e["file"] for e in results["errors"]]
        status["error_files"] = error_files[:5]  # Limit to 5 files
        status["has_errors"] = True
    else:
        status["has_errors"] = False
    
    if results["warnings"]:
        warning_files = [w["file"] for w in results["warnings"]]
        status["warning_files"] = warning_files[:5]  # Limit to 5 files
        status["has_warnings"] = True
    else:
        status["has_warnings"] = False
    
    return status

def main():
    try:
        # Get parameters from Snakemake
        log_files = [str(f) for f in snakemake.input.logs]
        error_patterns = snakemake.params.error_patterns
        warning_patterns = snakemake.params.warning_patterns
        output_report = snakemake.output.report
        output_status = snakemake.output.status
        
        # Log what we're doing
        log_with_context(logger, "INFO", f"Checking {len(log_files)} log files for errors and warnings")
        
        # Check logs
        results = check_logs_for_patterns(log_files, error_patterns, warning_patterns)
        
        # Generate status summary
        status = generate_status_summary(results)
        
        # Write results to files
        with open(output_report, 'w') as f:
            json.dump(results, f, indent=2)
        
        with open(output_status, 'w') as f:
            json.dump(status, f, indent=2)
        
        # Log summary
        if results["total_errors"] > 0:
            log_with_context(logger, "ERROR", f"Found {results['total_errors']} errors in {len(results['errors'])} log files")
        elif results["total_warnings"] > 0:
            log_with_context(logger, "WARNING", f"Found {results['total_warnings']} warnings in {len(results['warnings'])} log files")
        else:
            log_with_context(logger, "INFO", "No errors or warnings found in log files")
        
        return 0 if results["status"] == "ok" else 1
    
    except Exception as e:
        log_with_context(logger, "ERROR", f"Error checking logs: {e}")
        return 1

if __name__ == "__main__":
    if "snakemake" not in globals():
        logger.error("This script is designed to be run through Snakemake")
        sys.exit(1)
    
    sys.exit(main())