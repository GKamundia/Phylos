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
    start_time = time.time()
    # filepath: c:\Users\Anarchy\Documents\Data_Science\NextStrain\rvf-nextstrain\scripts\monitoring\monitor_duration.py
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
    start_time = time.time()
    