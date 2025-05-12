#!/usr/bin/env python3
"""
Send notifications about workflow completion via email or Slack
"""

import os
import json
import argparse
import logging
import smtplib
import requests
import yaml
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Send notifications about workflow completion")
    parser.add_argument("--log-file", help="Path to log file to analyze")
    parser.add_argument("--config", default="config/master_config.yaml", help="Path to master config file")
    return parser.parse_args()

def read_config(config_path):
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
            
        notification_config = config.get("common", {}).get("workflow", {}).get("notification", {})
        return {
            "enabled": notification_config.get("enabled", False),
            "email": notification_config.get("email", ""),
            "slack_webhook": notification_config.get("slack_webhook", ""),
            "pathogen": config.get("active_pathogen", "unknown"),
            "pathogen_name": config.get("pathogens", {}).get(config.get("active_pathogen", ""), {}).get("name", "Unknown Pathogen")
        }
    except Exception as e:
        logger.error(f"Failed to read config: {e}")
        return {
            "enabled": False,
            "email": "",
            "slack_webhook": "",
            "pathogen": "unknown",
            "pathogen_name": "Unknown Pathogen"
        }

def analyze_log(log_path):
    """Extract build status and key information from log file"""
    if not os.path.exists(log_path):
        return {"status": "unknown", "message": f"Log file not found: {log_path}"}
    
    try:
        with open(log_path, 'r') as f:
            log_content = f.read()
        
        # Check for success/error indicators
        if "ERROR:" in log_content:
            status = "failed"
            # Extract error message
            error_lines = [line for line in log_content.split('\n') if "ERROR:" in line]
            message = error_lines[-1] if error_lines else "Unknown error occurred"
        elif "build completed successfully" in log_content:
            status = "success"
            message = "Build completed successfully"
        else:
            status = "unknown"
            message = "Build status unclear"
        
        return {
            "status": status,
            "message": message,
            "log_path": log_path,
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
    except Exception as e:
        logger.error(f"Error analyzing log: {e}")
        return {"status": "error", "message": f"Failed to analyze log: {e}"}

def send_email_notification(recipient, subject, content):
    """Send notification email"""
    if not recipient:
        logger.warning("No email recipient specified, skipping email notification")
        return False
    
    try:
        # Create email
        msg = MIMEMultipart()
        msg["From"] = "nextstrain-notifications@example.com"  # Replace with actual email
        msg["To"] = recipient
        msg["Subject"] = subject
        
        # Attach message
        msg.attach(MIMEText(content, "plain"))
        
        # Connect to SMTP server
        # NOTE: In production, use environment variables or secure storage for credentials
        server = smtplib.SMTP("smtp.example.com", 587)  # Replace with actual SMTP server
        server.starttls()
        server.login("username", "password")  # Replace with actual credentials
        server.send_message(msg)
        server.quit()
        
        logger.info(f"Email notification sent to {recipient}")
        return True
    except Exception as e:
        logger.error(f"Failed to send email: {e}")
        return False

def send_slack_notification(webhook_url, content):
    """Send notification to Slack"""
    if not webhook_url:
        logger.warning("No Slack webhook URL provided, skipping Slack notification")
        return False
    
    try:
        # Create Slack message payload
        payload = {
            "text": content
        }
        
        # Send to Slack
        response = requests.post(
            webhook_url,
            json=payload,
            headers={"Content-Type": "application/json"}
        )
        
        if response.status_code == 200:
            logger.info("Slack notification sent successfully")
            return True
        else:
            logger.error(f"Failed to send Slack notification: HTTP {response.status_code}")
            return False
    except Exception as e:
        logger.error(f"Failed to send Slack notification: {e}")
        return False

def main():
    args = parse_args()
    config = read_config(args.config)
    
    if not config["enabled"]:
        logger.info("Notifications are disabled in config, exiting")
        return 0
    
    # Analyze build log if provided
    build_info = analyze_log(args.log_file) if args.log_file else {
        "status": "unknown", 
        "message": "No log file provided", 
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    
    # Format message
    status_emoji = "✅" if build_info["status"] == "success" else "❌" if build_info["status"] == "failed" else "⚠️"
    
    message = f"""
{status_emoji} Nextstrain Build Notification: {build_info["status"].upper()}

Pathogen: {config["pathogen_name"]} ({config["pathogen"]})
Status: {build_info["status"]}
Time: {build_info["timestamp"]}
Message: {build_info["message"]}

{f'Log file: {build_info["log_path"]}' if "log_path" in build_info else ''}
"""
    
    # Send notifications
    if config["email"]:
        subject = f"Nextstrain Build {build_info['status'].upper()}: {config['pathogen_name']}"
        send_email_notification(config["email"], subject, message)
    
    if config["slack_webhook"]:
        send_slack_notification(config["slack_webhook"], message)
    
    return 0

if __name__ == "__main__":
    main()