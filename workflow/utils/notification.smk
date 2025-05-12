"""
Rules for sending notifications about workflow progress and errors
"""

# Generate workflow summary
rule workflow_summary:
    input:
        # Use dynamic output from the all rule
        get_final_outputs
    output:
        summary = f"logs/summary/{output_prefix}_workflow_summary.json"
    params:
        pathogen = config["pathogen"],
        pathogen_name = config["pathogen_name"],
        segment_mode = segment_mode,
        segments = segments if segment_mode == "multi" else [config["data"].get("segment", "")]
    log:
        f"logs/workflow_summary_{output_prefix}.log"
    benchmark:
        f"benchmarks/workflow_summary_{output_prefix}.txt"
    resources:
        mem_mb = 1000
    script:
        "../../scripts/generate_workflow_summary.py"

# Send email notification (will be skipped if email settings are not configured)
rule send_email_notification:
    input:
        summary = f"logs/summary/{output_prefix}_workflow_summary.json"
    output:
        notification = f"logs/notifications/{output_prefix}_email_sent.txt"
    params:
        enabled = config.get("workflow", {}).get("notification", {}).get("enabled", False),
        email = config.get("workflow", {}).get("notification", {}).get("email", ""),
        subject = f"Nextstrain {config['pathogen_name']} Build Completed"
    log:
        f"logs/send_email_notification_{output_prefix}.log"
    resources:
        mem_mb = 1000
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.notification})
        
        if [ "{params.enabled}" = "True" ] && [ -n "{params.email}" ]; then
            python scripts/send_notification.py \
                --type email \
                --recipient "{params.email}" \
                --subject "{params.subject}" \
                --input-summary {input.summary} \
                --output {output.notification} \
                > {log} 2>&1
        else
            echo "Email notifications disabled or no email configured." > {output.notification}
            echo "Email notifications disabled or no email configured." > {log}
        fi
        """

# Send Slack notification (will be skipped if Slack settings are not configured)
rule send_slack_notification:
    input:
        summary = f"logs/summary/{output_prefix}_workflow_summary.json"
    output:
        notification = f"logs/notifications/{output_prefix}_slack_sent.txt"
    params:
        enabled = config.get("workflow", {}).get("notification", {}).get("enabled", False),
        webhook = config.get("workflow", {}).get("notification", {}).get("slack_webhook", "")
    log:
        f"logs/send_slack_notification_{output_prefix}.log"
    resources:
        mem_mb = 1000
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.notification})
        
        if [ "{params.enabled}" = "True" ] && [ -n "{params.webhook}" ]; then
            python scripts/send_notification.py \
                --type slack \
                --webhook "{params.webhook}" \
                --input-summary {input.summary} \
                --output {output.notification} \
                > {log} 2>&1
        else
            echo "Slack notifications disabled or no webhook configured." > {output.notification}
            echo "Slack notifications disabled or no webhook configured." > {log}
        fi
        """

# Combined notification rule for convenience
rule notify:
    input:
        email = f"logs/notifications/{output_prefix}_email_sent.txt",
        slack = f"logs/notifications/{output_prefix}_slack_sent.txt"
    output:
        touch(f"logs/notifications/{output_prefix}_all_sent.txt")