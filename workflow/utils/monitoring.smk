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

# Run alert system to check for issues and send notifications
rule run_alert_system:
    input:
        error_report = f"logs/monitoring/{output_prefix}_error_report.json",
        status = f"logs/monitoring/{output_prefix}_status.json",
        duration = f"logs/monitoring/{output_prefix}_duration_metrics.json"
    output:
        alerts = f"logs/monitoring/{output_prefix}_alerts.json"
    params:
        config_path = "config/master_config.yaml",
        logs_dir = "logs"
    log:
        f"logs/alert_system_{output_prefix}.log"
    resources:
        mem_mb = 1000
    shell:
        """
        python scripts/monitoring/alert_system.py \
            --config {params.config_path} \
            --logs-dir {params.logs_dir} \
            --output {output.alerts} \
            > {log} 2>&1
        """

# Generate comprehensive monitoring dashboard
rule generate_monitoring_dashboard:
    input:
        performance = f"logs/performance/{output_prefix}_performance_report.json",
        duration = f"logs/monitoring/{output_prefix}_duration_metrics.json",
        errors = f"logs/monitoring/{output_prefix}_error_report.json",
        alerts = f"logs/monitoring/{output_prefix}_alerts.json",
        summary = f"logs/summary/{output_prefix}_workflow_summary.json"
    output:
        dashboard = f"logs/monitoring/{output_prefix}_dashboard.html"
    params:
        pathogen = config["pathogen"],
        pathogen_name = config["pathogen_name"],
        segment_mode = segment_mode
    log:
        f"logs/monitoring_dashboard_{output_prefix}.log"
    resources:
        mem_mb = 1500
    script:
        "../../scripts/monitoring/generate_dashboard.py"