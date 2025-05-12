"""
Backup rules for Snakemake workflow
"""

# Create a backup of critical data
rule backup_critical_data:
    output:
        marker = "logs/backup/{date}_{prefix}_backup_complete.txt"
    params:
        backup_name = lambda wildcards: f"{wildcards.prefix}_{wildcards.date}"
    log:
        "logs/backup/backup_{date}_{prefix}.log"
    shell:
        """
        python scripts/backup_data.py backup \
            --name {params.backup_name} \
            > {log} 2>&1 && \
        echo "Backup completed at $(date)" > {output.marker}
        """

# Monitor backup health
rule check_backup_health:
    output:
        report = "logs/monitoring/{prefix}_backup_health.json"
    params:
        max_age = config.get("backup", {}).get("max_age", 7),
        email = config.get("backup", {}).get("alert_email", ""),
        backup_dir = config.get("backup", {}).get("directory", "backups")
    log:
        "logs/monitoring/backup_health_{prefix}.log"
    shell:
        """
        python scripts/monitoring/backup_monitor.py \
            --backup-dir {params.backup_dir} \
            --max-age {params.max_age} \
            --output {output.report} \
            --log-file {log} \
            {params.email and '--email ' + params.email or ''} \
            > {log} 2>&1
        """