#!/usr/bin/env python3
"""
Monitor workflow duration and alert if thresholds are exceeded
"""

import os
import sys
import json
import time
from datetime import datetime
from pathlib import Path

import snakemake

# Add parent directory to path for importing utils
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.log_utils import setup_logger, log_execution_stats, log_with_context

# Configure logger
logger = setup_logger(
    name="monitor_duration", 
    log_file=snakemake.log[0] if "snakemake" in globals() else "monitor_duration.log",
    level="INFO"
)

def parse_duration_from_summary(summary_file):
    """Parse duration information from workflow summary JSON"""
    try:
        with open(summary_file, 'r') as f:
            summary = json.load(f)
        
        runtime_seconds = summary.get("runtime", {}).get("seconds", 0)
        runtime_formatted = summary.get("runtime", {}).get("formatted", "unknown")
        
        return {
            "runtime_seconds": runtime_seconds,
            "runtime_formatted": runtime_formatted
        }
    
    except Exception as e:
        logger.error(f"Error parsing duration from summary file {summary_file}: {e}")
        return {
            "runtime_seconds": 0,
            "runtime_formatted": "unknown",
            "error": str(e)
        }

def check_duration_thresholds(duration, thresholds):
    """Check if duration exceeds warning or critical thresholds"""
    runtime_minutes = duration["runtime_seconds"] / 60
    
    status = "normal"
    message = ""
    
    if runtime_minutes > thresholds["critical"]:
        status = "critical"
        message = f"Runtime exceeds critical threshold of {thresholds['critical']} minutes"
    elif runtime_minutes > thresholds["warning"]:
        status = "warning"
        message = f"Runtime exceeds warning threshold of {thresholds['warning']} minutes"
    
    return {
        "status": status,
        "message": message,
        "runtime_minutes": runtime_minutes
    }

def generate_duration_report(duration, threshold_check, output_file, params):
    """Generate duration metrics report"""
    report = {
        "generated_at": datetime.now().isoformat(),
        "pathogen": params["pathogen"],
        "pathogen_name": params["pathogen_name"],
        "duration": {
            "seconds": duration["runtime_seconds"],
            "minutes": duration["runtime_seconds"] / 60,
            "formatted": duration["runtime_formatted"]
        },
        "thresholds": params["thresholds"],
        "status": threshold_check["status"],
        "message": threshold_check["message"]
    }
    
    # Write report
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    return report

def main():
    """
    Main function to monitor workflow durations and generate reports
    Designed to be run through Snakemake
    """
    start_time = time.time()
    
    try:
        # Get parameters from Snakemake
        summary_file = snakemake.input.summary
        output_file = snakemake.output.report
        params = {
            "pathogen": snakemake.params.pathogen,
            "pathogen_name": snakemake.params.pathogen_name,
            "thresholds": snakemake.params.thresholds
        }
        
        # Log the monitoring task start
        log_with_context(logger, "INFO", f"Monitoring workflow duration for {params['pathogen_name']}", {
            "summary_file": summary_file,
            "thresholds": params["thresholds"]
        })
        
        # Parse duration from summary file
        duration = parse_duration_from_summary(summary_file)
        if "error" in duration:
            log_with_context(logger, "ERROR", "Failed to parse workflow duration", {
                "summary_file": summary_file,
                "error": duration["error"]
            })
        else:
            # Log the parsed duration
            log_with_context(logger, "INFO", f"Workflow ran for {duration['runtime_formatted']}", {
                "seconds": duration["runtime_seconds"]
            })
        
        # Check against thresholds
        threshold_check = check_duration_thresholds(duration, params["thresholds"])
        
        # Log threshold check results
        if threshold_check["status"] != "normal":
            log_level = "WARNING" if threshold_check["status"] == "warning" else "ERROR"
            log_with_context(logger, log_level, threshold_check["message"], {
                "runtime_minutes": threshold_check["runtime_minutes"],
                "threshold": params["thresholds"][threshold_check["status"]]
            })
        else:
            log_with_context(logger, "INFO", "Workflow duration within normal thresholds", {
                "runtime_minutes": threshold_check["runtime_minutes"]
            })
        
        # Generate duration report
        report = generate_duration_report(duration, threshold_check, output_file, params)
        log_with_context(logger, "INFO", f"Generated duration report at {output_file}")
        
        # Log execution stats
        log_execution_stats(logger, start_time, {
            "operation": "monitor_duration",
            "pathogen": params["pathogen"],
            "status": threshold_check["status"],
            "runtime_seconds": duration["runtime_seconds"]
        })
        
        # Return exit code based on threshold check
        # This can be used by calling scripts to determine if further action is needed
        return 0 if threshold_check["status"] == "normal" else 1
    
    except Exception as e:
        log_with_context(logger, "CRITICAL", f"Error monitoring workflow duration: {str(e)}")
        log_execution_stats(logger, start_time, {
            "operation": "monitor_duration",
            "error": str(e)
        }, status="failed")
        return 1

if __name__ == "__main__":
    if "snakemake" not in globals():
        logger.error("This script is designed to be run through Snakemake")
        sys.exit(1)
    
    exit_code = main()
    sys.exit(exit_code)
