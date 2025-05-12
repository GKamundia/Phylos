#!/usr/bin/env python3
"""
Alert system for RVF-Nextstrain pipeline
Monitors logs and workflow status to send notifications for failures or warnings
"""

import os
import sys
import json
import re
import smtplib
import requests
import argparse
import logging
from datetime import datetime, timedelta
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Dict, List, Tuple, Optional, Any

# Add parent directory to path for importing utils
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.log_utils import setup_logger, log_with_context

# Configure logger
logger = setup_logger(name="alert_system", log_file="logs/monitoring/alert_system.log")

class AlertCondition:
    """Defines conditions that trigger alerts"""
    
    def __init__(
        self, 
        name: str,
        severity: str = "warning",
        pattern: Optional[str] = None,
        threshold: Optional[float] = None,
        window_minutes: int = 60
    ):
        self.name = name
        self.severity = severity  # warning, error, critical
        self.pattern = pattern    # regex pattern to match in logs
        self.threshold = threshold  # numeric threshold for metrics
        self.window_minutes = window_minutes  # time window to consider
    
    def to_dict(self) -> Dict:
        return {
            "name": self.name,
            "severity": self.severity,
            "pattern": self.pattern,
            "threshold": self.threshold,
            "window_minutes": self.window_minutes
        }

class AlertSystem:
    """System to detect issues and send notifications"""
    
    def __init__(
        self,
        config_path: str = "config/master_config.yaml",
        log_dir: str = "logs",
        default_conditions: Optional[List[AlertCondition]] = None
    ):
        self.config = self._load_config(config_path)
        self.log_dir = log_dir
        self.conditions = default_conditions or self._get_default_conditions()
    
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration from YAML file"""
        import yaml
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
            
            # Extract notification configuration
            notification_config = config.get("common", {}).get("workflow", {}).get("notification", {})
            return {
                "enabled": notification_config.get("enabled", False),
                "email": notification_config.get("email", ""),
                "slack_webhook": notification_config.get("slack_webhook", ""),
                "pathogen": config.get("active_pathogen", "unknown"),
                "pathogen_name": config.get("pathogens", {}).get(
                    config.get("active_pathogen", ""), {}
                ).get("name", "Unknown Pathogen")
            }
        except Exception as e:
            logger.error(f"Failed to load config: {e}")
            return {
                "enabled": False,
                "email": "",
                "slack_webhook": "",
                "pathogen": "unknown",
                "pathogen_name": "Unknown Pathogen"
            }
    
    def _get_default_conditions(self) -> List[AlertCondition]:
        """Define default alert conditions"""
        return [
            AlertCondition(
                name="critical_error",
                severity="critical",
                pattern=r"(CRITICAL|ERROR|Exception|failed|error:)",
                window_minutes=60
            ),
            AlertCondition(
                name="data_warning",
                severity="warning",
                pattern=r"(WARNING|missing|skipped|no sequences found)",
                window_minutes=60
            ),
            AlertCondition(
                name="long_runtime",
                severity="warning",
                threshold=120,  # Minutes
                window_minutes=240
            ),
            AlertCondition(
                name="failed_step",
                severity="error",
                pattern=r"(Error executing rule|MissingOutputException|WildcardError)",
                window_minutes=60
            )
        ]
    
    def scan_logs(self) -> List[Dict]:
        """Scan logs for alert conditions and return triggered alerts"""
        triggered_alerts = []
        
        # Set time window
        now = datetime.now()
        
        # Get all log files from the log directory
        log_files = []
        for root, _, files in os.walk(self.log_dir):
            for file in files:
                if file.endswith(".log"):
                    log_files.append(os.path.join(root, file))
        
        # Check each log file for pattern matches
        for log_file in log_files:
            try:
                # Skip if file is too old (based on modified time)
                file_mod_time = datetime.fromtimestamp(os.path.getmtime(log_file))
                if now - file_mod_time > timedelta(hours=24):
                    continue
                
                # Read and check file content
                with open(log_file, 'r') as f:
                    content = f.read()
                
                # Check each condition with a pattern
                for condition in self.conditions:
                    if condition.pattern and re.search(condition.pattern, content, re.IGNORECASE):
                        # Extract matching lines for context
                        matching_lines = []
                        for line in content.splitlines():
                            if re.search(condition.pattern, line, re.IGNORECASE):
                                matching_lines.append(line)
                        
                        # Limit to last 5 matching lines
                        matching_lines = matching_lines[-5:]
                        
                        # Add alert
                        triggered_alerts.append({
                            "condition": condition.to_dict(),
                            "file": log_file,
                            "timestamp": datetime.now().isoformat(),
                            "matching_lines": matching_lines
                        })
            except Exception as e:
                logger.error(f"Error scanning log file {log_file}: {e}")
        
        return triggered_alerts
    
    def check_performance_metrics(self) -> List[Dict]:
        """Check performance metrics for threshold violations"""
        triggered_alerts = []
        
        # Check duration metrics if available
        try:
            duration_files = []
            for root, _, files in os.walk(os.path.join(self.log_dir, "monitoring")):
                for file in files:
                    if "duration_metrics" in file and file.endswith(".json"):
                        duration_files.append(os.path.join(root, file))
            
            # Sort by modification time (newest first)
            duration_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
            
            if duration_files:
                # Check the most recent file
                with open(duration_files[0], 'r') as f:
                    metrics = json.load(f)
                
                # Check runtime threshold
                for condition in self.conditions:
                    if condition.threshold and "duration" in metrics:
                        runtime_minutes = metrics["duration"].get("minutes", 0)
                        if runtime_minutes > condition.threshold:
                            triggered_alerts.append({
                                "condition": condition.to_dict(),
                                "file": duration_files[0],
                                "timestamp": datetime.now().isoformat(),
                                "metric": "runtime",
                                "value": runtime_minutes,
                                "threshold": condition.threshold
                            })
        except Exception as e:
            logger.error(f"Error checking performance metrics: {e}")
        
        return triggered_alerts
    
    def send_email_notification(self, alerts: List[Dict]) -> bool:
        """Send email notification for alerts"""
        if not self.config["email"]:
            logger.warning("No email recipient configured, skipping email notification")
            return False
        
        try:
            # Prepare email content
            subject = f"[ALERT] Nextstrain {self.config['pathogen_name']} Pipeline Alert"
            
            # Count alerts by severity
            severity_count = {}
            for alert in alerts:
                severity = alert["condition"]["severity"]
                severity_count[severity] = severity_count.get(severity, 0) + 1
            
            # Create email body
            body = f"""
Nextstrain Pipeline Alert Report
===============================

Pathogen: {self.config['pathogen_name']} ({self.config['pathogen']})
Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Alert Summary:
{', '.join([f"{count} {sev}" for sev, count in severity_count.items()])}

Details:
"""
            
            # Add alert details
            for i, alert in enumerate(alerts, 1):
                body += f"\n{i}. [{alert['condition']['severity'].upper()}] {alert['condition']['name']}\n"
                body += f"   File: {alert['file']}\n"
                
                if "metric" in alert:
                    body += f"   Metric: {alert['metric']} = {alert['value']} (threshold: {alert['threshold']})\n"
                
                if "matching_lines" in alert and alert["matching_lines"]:
                    body += "   Context:\n"
                    for line in alert["matching_lines"]:
                        body += f"     {line.strip()}\n"
                
                body += "\n"
            
            body += """
Action Required:
Please check the logs and monitoring dashboards for more information.
"""
            
            # Create email
            msg = MIMEMultipart()
            msg["From"] = "nextstrain-alerts@example.com"
            msg["To"] = self.config["email"]
            msg["Subject"] = subject
            msg.attach(MIMEText(body, "plain"))
            
            # Configure SMTP details (in a real implementation, use env vars or secure config)
            smtp_server = "smtp.example.com"
            smtp_port = 587
            smtp_user = "username"
            smtp_pass = "password"
            
            # Send email
            server = smtplib.SMTP(smtp_server, smtp_port)
            server.starttls()
            server.login(smtp_user, smtp_pass)
            server.send_message(msg)
            server.quit()
            
            logger.info(f"Email alert sent to {self.config['email']}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to send email notification: {e}")
            return False
    
    def send_slack_notification(self, alerts: List[Dict]) -> bool:
        """Send Slack notification for alerts"""
        if not self.config["slack_webhook"]:
            logger.warning("No Slack webhook configured, skipping Slack notification")
            return False
        
        try:
            # Count alerts by severity
            severity_count = {}
            for alert in alerts:
                severity = alert["condition"]["severity"]
                severity_count[severity] = severity_count.get(severity, 0) + 1
            
            # Set emoji based on worst severity
            emoji = "⚠️"
            if "critical" in severity_count:
                emoji = "🚨"
            elif "error" in severity_count:
                emoji = "❌"
            
            # Create message blocks
            blocks = [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": f"{emoji} Nextstrain Pipeline Alert"
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": f"*Pathogen:* {self.config['pathogen_name']} ({self.config['pathogen']})\n*Time:* {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": f"*Alert Summary:*\n{', '.join([f'{count} {sev}' for sev, count in severity_count.items()])}"
                    }
                }
            ]
            
            # Add alert details (up to 5 to avoid too long message)
            for i, alert in enumerate(alerts[:5], 1):
                severity_emoji = "⚠️"
                if alert["condition"]["severity"] == "critical":
                    severity_emoji = "🚨"
                elif alert["condition"]["severity"] == "error":
                    severity_emoji = "❌"
                
                detail_text = f"{severity_emoji} *{alert['condition']['name']}*\n"
                detail_text += f"File: `{os.path.basename(alert['file'])}`\n"
                
                if "metric" in alert:
                    detail_text += f"Metric: {alert['metric']} = {alert['value']} (threshold: {alert['threshold']})\n"
                
                if "matching_lines" in alert and alert["matching_lines"]:
                    detail_text += "Context:\n"
                    for line in alert["matching_lines"][-2:]:  # Show just last 2 lines
                        detail_text += f"```{line.strip()}```\n"
                
                blocks.append({
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": detail_text
                    }
                })
            
            # Add separator
            blocks.append({"type": "divider"})
            
            # Add action footer
            blocks.append({
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": "Please check the logs and monitoring dashboards for more information."
                }
            })
            
            # Prepare payload
            payload = {
                "blocks": blocks
            }
            
            # Send to Slack
            response = requests.post(
                self.config["slack_webhook"],
                json=payload,
                headers={"Content-Type": "application/json"}
            )
            
            if response.status_code != 200:
                logger.error(f"Failed to send Slack notification: {response.status_code} {response.text}")
                return False
            
            logger.info("Slack notification sent successfully")
            return True
            
        except Exception as e:
            logger.error(f"Failed to send Slack notification: {e}")
            return False
    
    def run(self) -> Dict:
        """Run the alert system and send notifications if needed"""
        results = {
            "alerts_found": False,
            "alert_count": 0,
            "notifications_sent": False,
            "alerts": []
        }
        
        # Skip if notifications are disabled
        if not self.config["enabled"]:
            logger.info("Notifications are disabled in config")
            return results
        
        # Collect alerts
        log_alerts = self.scan_logs()
        metric_alerts = self.check_performance_metrics()
        all_alerts = log_alerts + metric_alerts
        
        if not all_alerts:
            logger.info("No alerts detected")
            return results
        
        # Update results
        results["alerts_found"] = True
        results["alert_count"] = len(all_alerts)
        results["alerts"] = all_alerts
        
        # Send notifications
        email_sent = False
        slack_sent = False
        
        if self.config["email"]:
            email_sent = self.send_email_notification(all_alerts)
        
        if self.config["slack_webhook"]:
            slack_sent = self.send_slack_notification(all_alerts)
        
        results["notifications_sent"] = email_sent or slack_sent
        
        return results

def parse_args():
    parser = argparse.ArgumentParser(description="Alert system for RVF-Nextstrain pipeline")
    parser.add_argument("--config", default="config/master_config.yaml", help="Path to config file")
    parser.add_argument("--logs-dir", default="logs", help="Directory containing log files")
    parser.add_argument("--output", help="Output file for alert report (JSON)")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Initialize alert system
    alert_system = AlertSystem(
        config_path=args.config,
        log_dir=args.logs_dir
    )
    
    # Run alert system
    results = alert_system.run()
    
    # Output results if requested
    if args.output:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
    
    # Return appropriate exit code
    if results["alerts_found"]:
        logger.warning(f"Found {results['alert_count']} alerts")
        return 1
    else:
        logger.info("No alerts detected")
        return 0

if __name__ == "__main__":
    sys.exit(main())