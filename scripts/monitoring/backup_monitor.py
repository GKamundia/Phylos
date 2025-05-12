
"""
Monitor backup health and send alerts for backup issues
"""

import os
import sys
import json
import datetime
import logging
import smtplib
import argparse
from email.mime.text import MIMEText
from pathlib import Path
from datetime import timedelta

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("backup_monitor")

def check_backup_health(backup_dir="backups", max_age_days=7):
    """
    Check the health of backups in the specified directory
    
    Args:
        backup_dir: Directory where backups are stored
        max_age_days: Maximum acceptable age of most recent backup in days
        
    Returns:
        dict: Health check results
    """
    backup_path = Path(backup_dir)
    
    if not backup_path.exists():
        return {
            "status": "error",
            "message": f"Backup directory {backup_dir} does not exist",
            "last_backup_age": None,
            "backup_count": 0,
            "issues": ["Backup directory not found"]
        }
    
    # Find all backups (directories with _backup_ in name or compressed archives)
    backups = []
    for item in backup_path.glob("*_backup_*"):
        if item.is_dir() or (item.is_file() and item.suffix in ['.gz', '.zip']):
            backups.append((item, item.stat().st_mtime))
    
    if not backups:
        return {
            "status": "error",
            "message": "No backups found",
            "last_backup_age": None,
            "backup_count": 0,
            "issues": ["No backups exist"]
        }
    
    # Sort by modification time (newest first)
    backups.sort(key=lambda x: x[1], reverse=True)
    newest_backup, newest_time = backups[0]
    
    # Calculate age of newest backup
    newest_date = datetime.datetime.fromtimestamp(newest_time)
    now = datetime.datetime.now()
    age_days = (now - newest_date).days
    
    issues = []
    status = "ok"
    message = "Backup health check passed"
    
    # Check if newest backup is too old
    if age_days > max_age_days:
        status = "warning"
        issues.append(f"Most recent backup is {age_days} days old (older than {max_age_days} days)")
        message = "Backup health check found issues"
    
    # Check if we have at least 2 backups
    if len(backups) < 2:
        if status != "error":  # Don't downgrade from error
            status = "warning"
        issues.append("Only one backup exists, redundancy is limited")
        message = "Backup health check found issues"
    
    # Check each backup's integrity
    for backup_path, _ in backups[:3]:  # Check only the 3 most recent backups
        if backup_path.is_dir():
            metadata_file = backup_path / "backup_metadata.json"
            if not metadata_file.exists():
                issues.append(f"Backup {backup_path.name} missing metadata file")
                status = "warning"
                message = "Backup health check found issues"
        elif backup_path.suffix == ".gz":
            # For compressed backups, we can't easily check integrity without extracting
            # Could implement more advanced checks like checking file size is reasonable
            pass
    
    return {
        "status": status,
        "message": message,
        "last_backup_age": age_days,
        "backup_count": len(backups),
        "last_backup": str(newest_backup),
        "last_backup_date": newest_date.isoformat(),
        "issues": issues
    }

def send_alert(health_check, email=None, log_file=None):
    """Send an alert if there are backup issues"""
    if health_check["status"] == "ok":
        return False
    
    alert_message = f"""
BACKUP MONITOR ALERT

Status: {health_check['status'].upper()}
Message: {health_check['message']}

Issues:
{chr(10).join('- ' + issue for issue in health_check['issues'])}

Last backup: {health_check.get('last_backup_date', 'Unknown')}
Age: {health_check.get('last_backup_age', 'Unknown')} days
Total backups: {health_check.get('backup_count', 0)}

Please check your backup system.
"""
    
    logger.warning(alert_message)
    
    # Write to log file if specified
    if log_file:
        with open(log_file, 'a') as f:
            f.write(f"\n{datetime.datetime.now().isoformat()}: BACKUP ALERT\n")
            f.write(alert_message)
            f.write("\n" + "-"*50 + "\n")
    
    # Send email if configured
    if email:
        try:
            msg = MIMEText(alert_message)
            msg['Subject'] = f"[RVF-Nextstrain] Backup {health_check['status'].upper()}: {health_check['message']}"
            msg['From'] = "nextstrain-monitor@example.com"
            msg['To'] = email
            
            s = smtplib.SMTP('localhost')
            s.send_message(msg)
            s.quit()
            logger.info(f"Alert email sent to {email}")
            return True
        except Exception as e:
            logger.error(f"Failed to send email alert: {e}")
    
    return True

def main():
    """Main function for the script"""
    parser = argparse.ArgumentParser(description="Monitor backup health")
    parser.add_argument("--backup-dir", default="backups", help="Directory containing backups")
    parser.add_argument("--max-age", type=int, default=7, help="Maximum acceptable age of most recent backup (days)")
    parser.add_argument("--email", help="Email address for alerts")
    parser.add_argument("--log-file", help="Log file for alerts")
    parser.add_argument("--output", help="Output file for health check results (JSON)")
    args = parser.parse_args()
    
    health_check = check_backup_health(args.backup_dir, args.max_age)
    
    # Print results
    print(f"Backup health status: {health_check['status']}")
    if health_check['issues']:
        print("Issues found:")
        for issue in health_check['issues']:
            print(f"  - {issue}")
    else:
        print("No issues found")
    print(f"Last backup: {health_check.get('last_backup_date', 'Unknown')}")
    print(f"Backup age: {health_check.get('last_backup_age', 'Unknown')} days")
    print(f"Total backups: {health_check['backup_count']}")
    
    # Save results to file if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(health_check, f, indent=2)
    
    # Send alert if issues were found
    if health_check["status"] != "ok":
        send_alert(health_check, args.email, args.log_file)
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())