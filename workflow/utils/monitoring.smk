"""
Monitoring and alerting rules for Nextstrain workflows
"""

# Generate performance report
rule generate_performance_report:
    input:
        benchmarks = lambda wildcards: list(Path("benchmarks").glob(f"*_{output_prefix}*"))
    output:
        report = f"logs/performance/{output_prefix}_performance_report.json",
        html = f"logs/performance/{output_prefix}_performance_dashboard.html"
    params:
        pathogen = config["pathogen"],
        pathogen_name = config["pathogen_name"],
        segment_mode = segment_mode,
        segments = segments if segment_mode == "multi" else [config["data"].get("segment", "")]
    log:
        f"logs/performance_report_{output_prefix}.log"
    resources:
        mem_mb = 1000
    script:
        "../../scripts/monitoring/generate_performance_report.py"

# Monitor workflow duration
rule monitor_workflow_duration:
    input:
        summary = f"logs/summary/{output_prefix}_workflow_summary.json"
    output:
        report = f"logs/monitoring/{output_prefix}_duration_metrics.json"
    params:
        pathogen = config["pathogen"],
        pathogen_name = config["pathogen_name"],
        thresholds = {
            "warning": 120,   # 2 hours in minutes
            "critical": 240   # 4 hours in minutes
        }
    log:
        f"logs/monitor_duration_{output_prefix}.log"
    resources:
        mem_mb = 1000
    script:
        "../../scripts/monitoring/monitor_duration.py"

# Check for execution errors
rule monitor_execution_errors:
    input:
        logs = lambda wildcards: list(Path("logs").glob("*.log"))
    output:
        report = f"logs/monitoring/{output_prefix}_error_report.json",
        status = f"logs/monitoring/{output_prefix}_status.json"
    params:
        error_patterns = [
            "ERROR", 
            "CRITICAL",
            "Exception",
            "failed",
            "error"
        ],
        warning_patterns = [
            "WARNING",
            "missing",
            "skipped",
            "Retry"
        ]
    log:
        f"logs/monitor_errors_{output_prefix}.log"
    resources:
        mem_mb = 1000
    script:
        "../../scripts/monitoring/check_logs_for_errors.py"

# Send monitoring alert
rule send_monitoring_alert:
    input:
        error_report = f"logs/monitoring/{output_prefix}_error_report.json",
        status = f"logs/monitoring/{output_prefix}_status.json",
        duration = f"logs/monitoring/{output_prefix}_duration_metrics.json"
    output:
        notification = f"logs/notifications/{output_prefix}_monitoring_alert.txt"
    params:
        enabled = config.get("workflow", {}).get("notification", {}).get("enabled", False),
        email = config.get("workflow", {}).get("notification", {}).get("email", ""),
        slack_webhook = config.get("workflow", {}).get("notification", {}).get("slack_webhook", ""),
        pathogen = config["pathogen"],
        pathogen_name = config["pathogen_name"]
    log:
        f"logs/send_monitoring_alert_{output_prefix}.log"
    resources:
        mem_mb = 1000
    script:
        "../../scripts/monitoring/send_monitoring_alert.py"