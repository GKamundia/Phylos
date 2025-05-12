#!/usr/bin/env python3
"""
Send notifications about workflow completion via email or Slack
"""

import os
import sys
import json
import argparse
import logging
import smtplib
import requests
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
    parser.add_argument("--type", required=True, choices=["email", "slack"], help="Notification type")
    parser.add_argument("--recipient", help="Email recipient")
    parser.add_argument("--subject", help="Email subject")
    parser.add_argument("--webhook", help="Slack webhook URL")
    parser.add_argument("--input-summary", required=True, help="Path to workflow summary JSON")
    parser.add_argument("--output", required=True, help="Output file to indicate success")
    return parser.parse_args()

def format_message(summary):
    """Format notification message based on workflow summary"""
    pathogen = summary.get("pathogen_name", "Unknown pathogen")
    timestamp = summary.get("timestamp", datetime.now().isoformat())
    segment_mode = summary.get("segment_mode", "unknown")
    segments = summary.get("segments", [])
    
    # Format sequence counts
    if segment_mode == "single":
        seq_counts = summary.get("sequence_counts", {})
        sequences_info = (
            f"Initial sequences: {seq_counts.get('raw', 0)}\n"
            f"Filtered: {seq_counts.get('filtered', 0)}\n"
            f"QC passed: {seq_counts.get('qc_passed', 0)}\n"
        )
    else:
        segment_counts = summary.get("segment_counts", {})
        sequences_info = f"Initial sequences: {summary.get('sequence_counts', {}).get('raw', 0)}\n"
        for segment, counts in segment_counts.items():
            sequences_info += (
                f"\nSegment {segment}:\n"
                f"  Filtered: {counts.get('filtered', 0)}\n"
                f"  QC passed: {counts.get('qc_passed', 0)}\n"
            )
    
    # Format runtime
    runtime = summary.get("runtime", {}).get("formatted", "unknown")
    
    # Create message
    message = (
        f"Nextstrain {pathogen} Build Completed\n"
        f"===============================\n\n"
        f"Timestamp: {timestamp}\n"
        f"Runtime: {runtime}\n"
        f"Segment mode: {segment_mode}\n"
        f"{'Segments: ' + ', '.join(segments) if segment_mode == 'multi' else ''}\n\n"
        f"Sequence Counts:\n"
        f"{sequences_info}\n"
    )
    
    return message

def send_email(recipient, subject, summary):
    """Send notification email"""
    try:
        message = format_message(summary)
        
        # Create email
        msg = MIMEMultipart()
        msg["From"] = "nextstrain-notifications@example.com"  # Replace with actual email
        msg["To"] = recipient
        msg["Subject"] = subject
        
        # Attach message
        msg.attach(MIMEText(message, "plain"))
        
        # Connect to SMTP server and send email
        # Note: In a real implementation, use environment variables or config
        # file for SMTP settings instead of hardcoding
        server = smtplib.SMTP("smtp.example.com", 587)  # Replace with actual SMTP server
        server.starttls()
        server.login("username", "password")  # Replace with actual credentials
        server.send_message(msg)
        server.quit()
        
        return True, "Email sent successfully"
        
    except Exception as e:
        return False, f"Failed to send email: {str(e)}"

def send_slack(webhook_url, summary):
    """Send notification to Slack"""
    try:
        message = format_message(summary)
        
        # Create Slack message payload
        payload = {
            "text": f"*Nextstrain Build Notification*\n```\n{message}\n```"
        }
        
        # Send to Slack
        response = requests.post(
            webhook_url,
            json=payload,
            headers={"Content-Type": "application/json"}
        )
        
        if response.status_code == 200:
            return True, "Slack notification sent successfully"
        else:
            return False, f"Failed to send Slack notification: HTTP {response.status_code}"
            
    except Exception as e:
        return False, f"Failed to send Slack notification: {str(e)}"

def main():
    args = parse_args()
    
    try:
        # Load workflow summary
        with open(args.input_summary, "r") as f:
            summary = json.load(f)
        
        # Send notification based on type
        if args.type == "email":
            if not args.recipient:
                logger.error("Email recipient is required for email notifications")
                return 1
                
            subject = args.subject or f"Nextstrain {summary.get('pathogen_name', 'Build')} Completed"
            success, message = send_email(args.recipient, subject, summary)
            
        elif args.type == "slack":
            if not args.webhook:
                logger.error("Webhook URL is required for Slack notifications")
                return 1
                
            success, message = send_slack(args.webhook, summary)
        
        else:
            logger.error(f"Unsupported notification type: {args.type}")
            return 1
        
        # Log result
        if success:
            logger.info(message)
            
            # Write output file to indicate success
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
            with open(args.output, "w") as f:
                f.write(f"{args.type} notification sent at {datetime.now().isoformat()}\n")
                
            return 0
        else:
            logger.error(message)
            return 1
            
    except Exception as e:
        logger.error(f"Error sending notification: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())