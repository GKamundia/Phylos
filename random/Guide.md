# 1.0 Comprehensive Implementation Plan for Reproducible Nextstrain Dashboards

This plan breaks down each phase of your #Guide.md into actionable steps.

Phase 1: Project Foundation and Environment Setup

1.1 Project Architecture Design
Review & Finalize Repository Structure:
Confirm the existing rvf-nextstrain directory structure (config, data, results, scripts, logs) aligns with the modular design principles in #Guide.md (Section 1.1).
Ensure clear separation of concerns between these directories.
Define Pathogen-Agnostic Workflow Principles:
Identify pathogen-specific parameters currently in config.yaml (e.g., search*term, segment, min_length).
Outline how these will be abstracted or made configurable for new pathogens (e.g., a master config pointing to pathogen-specific configs).
Establish Version Control Strategy:
Initialize a Git repository in the rvf-nextstrain root if not already done.
Review and update the .gitignore file (currently .gitignore) to include all transient files, large datasets (e.g., data/sequences/raw/*, results/\_), and environment-specific files, beyond just venv.
Decide on a strategy for Git LFS if any essential, non-text configuration files (e.g., reference genomes if large) need versioning.
1.2 Environment Configuration
Finalize Dockerfile for Reproducibility:
Create a Dockerfile in the rvf-nextstrain root.
Specify a base Nextstrain image (e.g., nextstrain/base:latest or a specific tagged version for better pinning).
Add RUN commands to install and pin versions of any additional Python packages or system dependencies required by your custom scripts (e.g., specific versions of biopython, pandas, pyyaml if different from the base image).
Document the chosen base image and pinned versions in your README.md.
Set Up Local Development Environment (Conda):
Create an environment.yml file in the rvf-nextstrain root for Conda.
List all necessary dependencies (Nextstrain CLI, Snakemake, Python, and specific versions of libraries used in your scripts) to mirror the Docker environment as closely as possible.
Document how to create and activate this Conda environment in the README.md.
Plan Cloud Infrastructure (Conceptual):
Identify potential future needs for cloud resources (S3 for data/results storage, EC2/Batch for computation) as outlined in #Guide.md (Section 1.2). This is for planning at this stage.
Phase 2: Data Acquisition and Curation System

2.1 Robust Data Integration
Enhance Data Fetching (download_sequences.py):
Review download_sequences.py. Plan modifications to make it more configurable for different NCBI search terms (beyond config.yaml) or even different data sources (e.g., by adding parameters or abstracting the fetching logic).
Consider how to handle different segments of RVF more explicitly if the current script doesn't fully address this based on config.yaml's segment setting.
Design Local/Remote Sequence Archive Strategy:
Decide whether to implement a versioned local archive (e.g., dated subdirectories within data/sequences/raw_archive/) or plan for S3 integration for raw downloaded data, as suggested in #Guide.md (Section 2.1).
Outline how scripts will interact with this archive.
Plan Incremental Data Fetching:
Conceptualize how download_sequences.py could be modified to fetch only new or updated sequences since the last run (e.g., by tracking last fetch dates or querying NCBI with date ranges).
2.2 Metadata Standardization Engine
Define and Implement Strict Metadata Schema:
Create a formal schema definition for your metadata (e.g., a JSON Schema file stored in config/metadata_schema.json). This schema should list all required and optional fields, their types, and any constraints.
Plan to integrate validation against this schema into prepare_metadata.py or a new validation script.
Enhance Metadata Preparation (prepare_metadata.py):
Review prepare_metadata.py. Ensure its date standardization is robust.
Plan additions for location standardization (e.g., using lat_longs.tsv more effectively, or integrating a geocoding library for missing coordinates if feasible).
Plan for host name normalization and other pathogen-specific field standardizations.
Design Metadata Enhancement Pipeline:
Outline steps for inferring missing critical metadata where possible (e.g., deriving division from location if a mapping exists).
Plan how to flag or report records with insufficient metadata for analysis.
2.3 Quality Control System
Design Multi-Level QC Pipeline:
Identify QC steps to be added to the Snakefile. This includes:
Initial filtering (already in Snakefile via augur filter, review parameters in config.yaml).
Plan for integrating Nextclade for sequence QC, clade assignment, and outlier detection. This would involve adding a new rule to the Snakefile and processing Nextclade's output.
Plan QC Reporting:
Determine what QC metrics to capture (e.g., number of sequences filtered at each step, reasons for exclusion).
Plan how these metrics will be logged or reported (e.g., summary files in results/qc_reports/).
Phase 3: Nextstrain Build Configuration and Optimization

3.1 Efficient Nextstrain Build Setup
Adapt Snakefile for Templating:
Review the existing Snakefile. Identify parts that are RVF-specific and could be parameterized further using config.yaml or additional pathogen-specific config files.
Ensure the Snakefile follows the Nextstrain pathogen-repo-guide structure as much as possible.
Implement Configurable Genome Segment Handling:
Verify that the segment parameter in config.yaml is effectively used by download_sequences.py and potentially by Snakefile rules if segment-specific analysis paths are needed (e.g., different reference sequences for alignment based on segment).
Plan how to handle "all" segments if that implies concatenation or separate analyses followed by a combined export.
Optimize Resource Usage in Snakefile:
Review threads parameters in Snakefile rules (align, tree).
Consider if other rules could benefit from thread/core specification.
Plan for smart caching of intermediate results (Snakemake does this by default, but ensure rules are defined to maximize this).
3.2 Advanced Snakemake Configuration
Refine Snakefile Modularity and Conditionality:
Break down complex Snakefile rules into smaller, more manageable ones if necessary.
Plan for conditional execution of rules based on config.yaml parameters (e.g., skipping a segment-specific rule if segment: "all" implies a different workflow).
Manage Dependencies:
Ensure the Dockerfile and environment.yml (Conda) accurately reflect all software dependencies and their pinned versions for Augur, Auspice, and any tools used by your scripts.
Implement Workflow Checkpointing:
Identify long-running or critical rules in the Snakefile where Snakemake checkpoints could be beneficial for resuming failed runs. (Snakemake's default behavior often handles this well if intermediate files are correctly defined).
Phase 4: Automation and Observability Framework

4.1 Reliable Scheduling System
Plan Workflow Scheduling (GitHub Actions / Cron):
Decide on the scheduling mechanism (GitHub Actions if the project will be on GitHub, or local cron jobs/Windows Task Scheduler for local automation).
Outline the script/command to be scheduled (e.g., nextstrain build . --configfile config/pathogen_X_config.yaml).
Design Update Coordination (if applicable):
If different parts of the pipeline (e.g., data download vs. full analysis) need to run on different schedules, plan how these will be coordinated.
4.2 Comprehensive Monitoring
Implement Structured Logging:
Review logging in download_sequences.py and prepare_metadata.py. Ensure they output informative messages.
Ensure all Snakefile rules redirect stdout and stderr to distinct log files in the logs/ directory (already implemented).
Plan Alerting System:
Conceptualize how failures in scheduled runs would trigger notifications (e.g., email via a cron job's mail capabilities, or using GitHub Actions notification features).
Plan Performance Tracking:
Consider adding simple timing commands (e.g., date before and after key shell commands in Snakefile) or using Snakemake's benchmarking capabilities to track rule runtimes.
4.3 Failure Recovery Systems
Design Retry Mechanisms:
For data download steps, consider adding simple retry logic within the Python script for transient network errors.
Snakemake can be re-run, and it will pick up from the last successful step.
Plan Manual Intervention Capabilities:
Ensure the Snakefile can be run to regenerate specific intermediate files or from specific points if needed (Snakemake's target-based execution).
Design Data Recovery Procedures:
Outline a strategy for backing up critical raw data and configurations. If using Git, regular pushes to a remote serve as a backup for code and small configs.
Phase 5: Dashboard and Visualization

5.1 Efficient Auspice Integration
Customize auspice_config.json for RVF:
Review and refine auspice_config.json. Ensure colorings, filters, and panels are optimal for RVF data.
Add pathogen-specific descriptions or annotations if needed.
Plan Auspice Performance Optimizations:
Consider the subsampling strategies defined in config.yaml and implemented via augur filter in the Snakefile as the primary way to manage dataset size for Auspice.
Design Metadata-Driven Auspice Configuration:
Conceptualize how auspice_config.json might be dynamically adjusted or selected based on the pathogen if you were to have multiple pathogen builds in the same overarching structure. For now, focus on making the RVF one robust.
5.2 Enhanced Dashboard Features (Conceptual)
Plan for Custom Visualization Components:
Identify if any RVF-specific visualizations are needed beyond standard Auspice panels. This is for future consideration.
Plan Server-Side Caching:
If deploying the dashboard to a high-traffic environment in the future, consider how server-side caching of the Auspice JSON might be implemented.
Plan Embedding Capabilities:
Consider if parts of the dashboard might need to be embedded elsewhere in the future.
Phase 6: Operationalization and Maintenance

6.1 Documentation and Knowledge Transfer
Update README.md Comprehensively:
Expand the existing README.md to include:
Detailed setup instructions for both Docker and local Conda environments.
Instructions on how to run the pipeline (e.g., nextstrain build .).
Explanation of the config.yaml parameters.
Description of the project structure.
Troubleshooting common issues.
Add Inline Documentation:
Ensure all Python scripts (scripts/\*.py) have clear docstrings and comments.
Add comments to the Snakefile explaining complex rules or logic.
Plan Onboarding Materials:
Outline what a new user/developer would need to know to understand and contribute to the project. The README.md will be the primary document for this.
6.2 Testing and Validation
Plan Automated Tests:
Unit Tests: Plan for simple unit tests for critical functions in your Python scripts (e.g., date parsing in prepare_metadata.py).
Integration Tests: Plan to use a small, static test dataset (sequences and metadata) to run the entire Snakefile workflow to ensure all steps integrate correctly. This test dataset could live in test_data/.
Smoke Tests: The integration test can also serve as a smoke test.
Set Up Continuous Integration (CI) (Future):
If the project moves to GitHub, plan to set up GitHub Actions to automatically run tests on pushes or pull requests.
Plan Data Validation Procedures:
Incorporate metadata schema validation (from step 2.2.1) into the workflow.
Plan for checks on the output Auspice JSON (e.g., ensuring it's valid JSON and contains expected top-level keys).
6.3 Extensibility for New Pathogens
Document Adaptation Process:
Create a new Markdown file (e.g., ADAPTING_FOR_NEW_PATHOGENS.md) outlining:
Which configuration files need to be created/modified (e.g., a new config/cholera_config.yaml, config/cholera_auspice_config.json, config/cholera_lat_longs.tsv).
Which scripts might need pathogen-specific logic (and how to add it, e.g., via conditional statements based on a config value).
How to create new Snakefile rules or modify existing ones if the core bioinformatics steps differ significantly.
Create Configuration Templates:
Use the existing RVF config files (config.yaml, auspice_config.json, lat_longs.tsv) as de facto templates. The adaptation guide will instruct users to copy and modify these.
Design Plugin System (Conceptual):
For significant variations in data processing or analysis for different pathogens, consider a more formal plugin-like structure for scripts or Snakefile rule sets in the future. For now, focus on clear separation and configuration.
Phase 7: Deployment and Launch

7.1 Production Deployment (Local Focus First)
Set Up "Staging" (Local Testing):
Your current local setup where you run nextstrain build . and nextstrain view results/auspice/ serves as your initial testing/staging environment.
Plan "Production" Environment Configuration (Local):
Define a clear, stable directory path for your automated runs if different from your development path.
Ensure the user account running scheduled tasks has necessary permissions.
Create Launch Plan (for RVF dashboard):
Steps to perform a full, clean run.
How to view the output.
How to verify the first automated update.
7.2 Ongoing Operations
Establish Update Procedures:
Document how to manually trigger an update if needed.
Document how to check logs for scheduled updates.
Plan Performance Optimization:
Periodically review runtimes (from logs or Snakemake reports) and dataset sizes to see if subsampling strategies in config.yaml need adjustment.
Implement User Feedback Loop (Conceptual):
If others will use this dashboard, think about how you would collect feedback.
This detailed plan should guide you through implementing your #Guide.md effectively, leveraging your existing project structure and focusing on the core goals of reproducibility and adaptability. Remember to tackle these steps iteratively.

# 2.0 Comprehensive Implementation Plan for Reproducible Nextstrain Dashboards

This plan breaks down each phase of your #Guide.md into actionable steps.

Phase 1: Project Foundation and Environment Setup

1.1 Project Architecture Design
Review & Finalize Repository Structure:
Confirm the existing rvf-nextstrain directory structure (config, data, results, scripts, logs) aligns with the modular design principles in #Guide.md (Section 1.1).
Ensure clear separation of concerns between these directories.
Define Pathogen-Agnostic Workflow Principles:
Identify pathogen-specific parameters currently in config.yaml (e.g., search*term, segment, min_length).
Outline how these will be abstracted or made configurable for new pathogens (e.g., a master config pointing to pathogen-specific configs).
Establish Version Control Strategy:
Initialize a Git repository in the rvf-nextstrain root if not already done.
Review and update the .gitignore file (currently .gitignore) to include all transient files, large datasets (e.g., data/sequences/raw/*, results/\_), and environment-specific files, beyond just venv.
Decide on a strategy for Git LFS if any essential, non-text configuration files (e.g., reference genomes if large) need versioning.
1.2 Environment Configuration
Finalize Dockerfile for Reproducibility:
Create a Dockerfile in the rvf-nextstrain root.
Specify a base Nextstrain image (e.g., nextstrain/base:latest or a specific tagged version for better pinning).
Add RUN commands to install and pin versions of any additional Python packages or system dependencies required by your custom scripts (e.g., specific versions of biopython, pandas, pyyaml if different from the base image).
Document the chosen base image and pinned versions in your README.md.
Set Up Local Development Environment (Conda):
Create an environment.yml file in the rvf-nextstrain root for Conda.
List all necessary dependencies (Nextstrain CLI, Snakemake, Python, and specific versions of libraries used in your scripts) to mirror the Docker environment as closely as possible.
Document how to create and activate this Conda environment in the README.md.
Plan Cloud Infrastructure (Conceptual):
Identify potential future needs for cloud resources (S3 for data/results storage, EC2/Batch for computation) as outlined in #Guide.md (Section 1.2). This is for planning at this stage.
Phase 2: Data Acquisition and Curation System

2.1 Robust Data Integration
Enhance Data Fetching (download_sequences.py):
Review download_sequences.py. Plan modifications to make it more configurable for different NCBI search terms (beyond config.yaml) or even different data sources (e.g., by adding parameters or abstracting the fetching logic).
Consider how to handle different segments of RVF more explicitly if the current script doesn't fully address this based on config.yaml's segment setting.
Design Local/Remote Sequence Archive Strategy:
Decide whether to implement a versioned local archive (e.g., dated subdirectories within data/sequences/raw_archive/) or plan for S3 integration for raw downloaded data, as suggested in #Guide.md (Section 2.1).
Outline how scripts will interact with this archive.
Plan Incremental Data Fetching:
Conceptualize how download_sequences.py could be modified to fetch only new or updated sequences since the last run (e.g., by tracking last fetch dates or querying NCBI with date ranges).
2.2 Metadata Standardization Engine
Define and Implement Strict Metadata Schema:
Create a formal schema definition for your metadata (e.g., a JSON Schema file stored in config/metadata_schema.json). This schema should list all required and optional fields, their types, and any constraints.
Plan to integrate validation against this schema into prepare_metadata.py or a new validation script.
Enhance Metadata Preparation (prepare_metadata.py):
Review prepare_metadata.py. Ensure its date standardization is robust.
Plan additions for location standardization (e.g., using lat_longs.tsv more effectively, or integrating a geocoding library for missing coordinates if feasible).
Plan for host name normalization and other pathogen-specific field standardizations.
Design Metadata Enhancement Pipeline:
Outline steps for inferring missing critical metadata where possible (e.g., deriving division from location if a mapping exists).
Plan how to flag or report records with insufficient metadata for analysis.
2.3 Quality Control System
Design Multi-Level QC Pipeline:
Identify QC steps to be added to the Snakefile. This includes:
Initial filtering (already in Snakefile via augur filter, review parameters in config.yaml).
Plan for integrating Nextclade for sequence QC, clade assignment, and outlier detection. This would involve adding a new rule to the Snakefile and processing Nextclade's output.
Plan QC Reporting:
Determine what QC metrics to capture (e.g., number of sequences filtered at each step, reasons for exclusion).
Plan how these metrics will be logged or reported (e.g., summary files in results/qc_reports/).
Phase 3: Nextstrain Build Configuration and Optimization

3.1 Efficient Nextstrain Build Setup
Adapt Snakefile for Templating:
Review the existing Snakefile. Identify parts that are RVF-specific and could be parameterized further using config.yaml or additional pathogen-specific config files.
Ensure the Snakefile follows the Nextstrain pathogen-repo-guide structure as much as possible.
Implement Configurable Genome Segment Handling:
Verify that the segment parameter in config.yaml is effectively used by download_sequences.py and potentially by Snakefile rules if segment-specific analysis paths are needed (e.g., different reference sequences for alignment based on segment).
Plan how to handle "all" segments if that implies concatenation or separate analyses followed by a combined export.
Optimize Resource Usage in Snakefile:
Review threads parameters in Snakefile rules (align, tree).
Consider if other rules could benefit from thread/core specification.
Plan for smart caching of intermediate results (Snakemake does this by default, but ensure rules are defined to maximize this).
3.2 Advanced Snakemake Configuration
Refine Snakefile Modularity and Conditionality:
Break down complex Snakefile rules into smaller, more manageable ones if necessary.
Plan for conditional execution of rules based on config.yaml parameters (e.g., skipping a segment-specific rule if segment: "all" implies a different workflow).
Manage Dependencies:
Ensure the Dockerfile and environment.yml (Conda) accurately reflect all software dependencies and their pinned versions for Augur, Auspice, and any tools used by your scripts.
Implement Workflow Checkpointing:
Identify long-running or critical rules in the Snakefile where Snakemake checkpoints could be beneficial for resuming failed runs. (Snakemake's default behavior often handles this well if intermediate files are correctly defined).
Phase 4: Automation and Observability Framework

4.1 Reliable Scheduling System
Plan Workflow Scheduling (GitHub Actions / Cron):
Decide on the scheduling mechanism (GitHub Actions if the project will be on GitHub, or local cron jobs/Windows Task Scheduler for local automation).
Outline the script/command to be scheduled (e.g., nextstrain build . --configfile config/pathogen_X_config.yaml).
Design Update Coordination (if applicable):
If different parts of the pipeline (e.g., data download vs. full analysis) need to run on different schedules, plan how these will be coordinated.
4.2 Comprehensive Monitoring
Implement Structured Logging:
Review logging in download_sequences.py and prepare_metadata.py. Ensure they output informative messages.
Ensure all Snakefile rules redirect stdout and stderr to distinct log files in the logs/ directory (already implemented).
Plan Alerting System:
Conceptualize how failures in scheduled runs would trigger notifications (e.g., email via a cron job's mail capabilities, or using GitHub Actions notification features).
Plan Performance Tracking:
Consider adding simple timing commands (e.g., date before and after key shell commands in Snakefile) or using Snakemake's benchmarking capabilities to track rule runtimes.
4.3 Failure Recovery Systems
Design Retry Mechanisms:
For data download steps, consider adding simple retry logic within the Python script for transient network errors.
Snakemake can be re-run, and it will pick up from the last successful step.
Plan Manual Intervention Capabilities:
Ensure the Snakefile can be run to regenerate specific intermediate files or from specific points if needed (Snakemake's target-based execution).
Design Data Recovery Procedures:
Outline a strategy for backing up critical raw data and configurations. If using Git, regular pushes to a remote serve as a backup for code and small configs.
Phase 5: Dashboard and Visualization

5.1 Efficient Auspice Integration
Customize auspice_config.json for RVF:
Review and refine auspice_config.json. Ensure colorings, filters, and panels are optimal for RVF data.
Add pathogen-specific descriptions or annotations if needed.
Plan Auspice Performance Optimizations:
Consider the subsampling strategies defined in config.yaml and implemented via augur filter in the Snakefile as the primary way to manage dataset size for Auspice.
Design Metadata-Driven Auspice Configuration:
Conceptualize how auspice_config.json might be dynamically adjusted or selected based on the pathogen if you were to have multiple pathogen builds in the same overarching structure. For now, focus on making the RVF one robust.
5.2 Enhanced Dashboard Features (Conceptual)
Plan for Custom Visualization Components:
Identify if any RVF-specific visualizations are needed beyond standard Auspice panels. This is for future consideration.
Plan Server-Side Caching:
If deploying the dashboard to a high-traffic environment in the future, consider how server-side caching of the Auspice JSON might be implemented.
Plan Embedding Capabilities:
Consider if parts of the dashboard might need to be embedded elsewhere in the future.
Phase 6: Operationalization and Maintenance

6.1 Documentation and Knowledge Transfer
Update README.md Comprehensively:
Expand the existing README.md to include:
Detailed setup instructions for both Docker and local Conda environments.
Instructions on how to run the pipeline (e.g., nextstrain build .).
Explanation of the config.yaml parameters.
Description of the project structure.
Troubleshooting common issues.
Add Inline Documentation:
Ensure all Python scripts (scripts/\*.py) have clear docstrings and comments.
Add comments to the Snakefile explaining complex rules or logic.
Plan Onboarding Materials:
Outline what a new user/developer would need to know to understand and contribute to the project. The README.md will be the primary document for this.
6.2 Testing and Validation
Plan Automated Tests:
Unit Tests: Plan for simple unit tests for critical functions in your Python scripts (e.g., date parsing in prepare_metadata.py).
Integration Tests: Plan to use a small, static test dataset (sequences and metadata) to run the entire Snakefile workflow to ensure all steps integrate correctly. This test dataset could live in test_data/.
Smoke Tests: The integration test can also serve as a smoke test.
Set Up Continuous Integration (CI) (Future):
If the project moves to GitHub, plan to set up GitHub Actions to automatically run tests on pushes or pull requests.
Plan Data Validation Procedures:
Incorporate metadata schema validation (from step 2.2.1) into the workflow.
Plan for checks on the output Auspice JSON (e.g., ensuring it's valid JSON and contains expected top-level keys).
6.3 Extensibility for New Pathogens
Document Adaptation Process:
Create a new Markdown file (e.g., ADAPTING_FOR_NEW_PATHOGENS.md) outlining:
Which configuration files need to be created/modified (e.g., a new config/cholera_config.yaml, config/cholera_auspice_config.json, config/cholera_lat_longs.tsv).
Which scripts might need pathogen-specific logic (and how to add it, e.g., via conditional statements based on a config value).
How to create new Snakefile rules or modify existing ones if the core bioinformatics steps differ significantly.
Create Configuration Templates:
Use the existing RVF config files (config.yaml, auspice_config.json, lat_longs.tsv) as de facto templates. The adaptation guide will instruct users to copy and modify these.
Design Plugin System (Conceptual):
For significant variations in data processing or analysis for different pathogens, consider a more formal plugin-like structure for scripts or Snakefile rule sets in the future. For now, focus on clear separation and configuration.
Phase 7: Deployment and Launch

7.1 Production Deployment (Local Focus First)
Set Up "Staging" (Local Testing):
Your current local setup where you run nextstrain build . and nextstrain view results/auspice/ serves as your initial testing/staging environment.
Plan "Production" Environment Configuration (Local):
Define a clear, stable directory path for your automated runs if different from your development path.
Ensure the user account running scheduled tasks has necessary permissions.
Create Launch Plan (for RVF dashboard):
Steps to perform a full, clean run.
How to view the output.
How to verify the first automated update.
7.2 Ongoing Operations
Establish Update Procedures:
Document how to manually trigger an update if needed.
Document how to check logs for scheduled updates.
Plan Performance Optimization:
Periodically review runtimes (from logs or Snakemake reports) and dataset sizes to see if subsampling strategies in config.yaml need adjustment.
Implement User Feedback Loop (Conceptual):
If others will use this dashboard, think about how you would collect feedback.
This detailed plan should guide you through implementing your #Guide.md effectively, leveraging your existing project structure and focusing on the core goals of reproducibility and adaptability. Remember to tackle these steps iteratively.
