# Setting Up the Scheduling System

This document explains how to set up automated scheduling for the Nextstrain RVF pipeline.

## GitHub Actions (Cloud-based Automation)

If your repository is hosted on GitHub:

1. Make sure your repository contains the `.github/workflows/scheduled_builds.yml` file
2. Push your changes to GitHub
3. The workflow will automatically run on schedule (default: every Monday at 2:00 AM UTC)
4. You can also manually trigger the workflow from the "Actions" tab in your GitHub repository

### Customizing the Schedule

To change the schedule, edit the `cron` expression in `.github/workflows/scheduled_builds.yml`:

```yaml
on:
  schedule:
    - cron: "0 2 * * 1" # Every Monday at 2:00 AM UTC
```

Windows Task Scheduler (Local Automation)
To set up automated runs on a Windows system:

Open Task Scheduler (search for it in the Start menu)
Click "Create Basic Task..."
Name: "Nextstrain RVF Weekly Build" (or whatever you prefer)
Trigger: Choose when you want the task to run (e.g., Weekly)
Action: "Start a program"
Program/script: Browse to scripts\run_scheduled_build.bat
Add arguments (optional):
For full build: leave empty
For data-only update: true
For specific pathogen: false cholera (replace "cholera" with any pathogen in your config)
Complete the wizard
For Different Update Schedules
If you need separate schedules for data acquisition vs. full analysis:

Create a task as above for data acquisition, with argument true
Create another task for full analysis (leave arguments empty)
Set different schedules for each task
