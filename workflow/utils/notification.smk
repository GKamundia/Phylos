"""
Rules and functions for notifications on workflow completion and failures
"""

import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

# Import or redefine get_final_outputs function to avoid the NameError
# We'll redefine it here to avoid circular imports
def get_notification_outputs(wildcards=None):
    """Generate a list of final outputs based on configuration for notifications"""
    outputs = []
    
    # QC reports always included
    outputs.append(f"results/qc_reports/{output_prefix}_qc_summary.json")
    
    # Single segment mode outputs
    if segment_mode == "single":
        outputs.append(f"results/auspice/{output_prefix}.json")
    
    # Multi-segment mode outputs
    else:
        # Individual segment outputs
        for segment in segments:
            outputs.append(f"results/segments/{segment}/auspice/{output_prefix}_{segment}.json")
        
        # Combined visualization if enabled
        if config.get("workflow", {}).get("create_combined_view", True):
            outputs.append(f"results/auspice/{output_prefix}_combined.json")
    
    return outputs

# Notification rule that runs after successful completion
rule send_notification:
    input:
        get_notification_outputs  # Use the local function instead of the Snakefile's function
    output:
        notification_log = "logs/notifications/notification_{timestamp}.log"
    params:
        pathogen = config["pathogen"],
        pathogen_name = config["pathogen_name"],
        email = config.get("workflow", {}).get("notification", {}).get("email", ""),
        slack_webhook = config.get("workflow", {}).get("notification", {}).get("slack_webhook", ""),
        enabled = config.get("workflow", {}).get("notification", {}).get("enabled", False)
    run:
        import time
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        
        if not params.enabled:
            with open(output.notification_log, "w") as f:
                f.write(f"Notifications disabled in config. Skipping.\n")
            print("Notifications disabled in config. Skipping.")
            return
        
        # Prepare notification message
        message = f"""
        Nextstrain Build Completed: {params.pathogen_name} ({params.pathogen})
        
        Build completed successfully on {time.strftime("%Y-%m-%d %H:%M:%S")}
        
        Generated outputs:
        """
        
        for file_path in input:
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path) / (1024 * 1024)  # Convert to MB
                message += f"- {file_path} ({file_size:.2f} MB)\n"
        
        # Log the notification
        with open(output.notification_log, "w") as f:
            f.write(f"Notification sent at {timestamp}\n")
            f.write(message)
        
        # Send email notification if configured
        if params.email:
            try:
                send_email_notification(params.email, f"Nextstrain Build Complete: {params.pathogen_name}", message)
                with open(output.notification_log, "a") as f:
                    f.write(f"\nEmail notification sent to {params.email}")
            except Exception as e:
                with open(output.notification_log, "a") as f:
                    f.write(f"\nFailed to send email notification: {str(e)}")
        
        # Send Slack notification if configured
        if params.slack_webhook:
            try:
                send_slack_notification(params.slack_webhook, message)
                with open(output.notification_log, "a") as f:
                    f.write(f"\nSlack notification sent via webhook")
            except Exception as e:
                with open(output.notification_log, "a") as f:
                    f.write(f"\nFailed to send Slack notification: {str(e)}")

# Helper functions for notifications
def send_email_notification(recipient, subject, message):
    """Send email notification"""
    # This is a placeholder - implementation would depend on your email setup
    # You might use smtplib for a simple solution or a dedicated email service
    pass

def send_slack_notification(webhook_url, message):
    """Send Slack notification via webhook"""
    # This is a placeholder - implementation would use requests to post to Slack webhook
    pass

# Configure onstart/onsuccess/onerror hooks if needed
onsuccess:
    if config.get("workflow", {}).get("notification", {}).get("enabled", False):
        print(f"Build completed successfully for {config['pathogen_name']} ({config['pathogen']})")
        # Notification is handled by the rule above when using Snakemake's standard execution

onerror:
    if config.get("workflow", {}).get("notification", {}).get("enabled", False):
        # For error notification, we can't rely on the rule since it won't run if there's an error
        print(f"Build failed for {config['pathogen_name']} ({config['pathogen']})")
        # Here you would add code to send failure notifications