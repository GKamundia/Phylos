# Rift Valley Fever Nextstrain Dashboard

This repository contains a complete workflow for building a Nextstrain-based genomic surveillance dashboard for Rift Valley Fever (RVF) virus. The pipeline is designed to be easily adaptable for other pathogens.

## Setup Instructions

### Prerequisites

- Docker installed
- Python 3.7+
- Nextstrain CLI

### Installation

1. **Clone this repository:**

   ```bash
   git clone https://github.com/yourusername/rvf-nextstrain.git
   cd rvf-nextstrain
   ```

2. **Download data from NCBI:**

```bash
python scripts/download_ncbi_rvf_data.py --email your.email@example.com --segment L --max-sequences 500
OR
python scripts/download_ncbi_rvf_data.py --email your.email@example.com --segment L --max-sequences 500 --api-key YOUR_API_KEY
```

3. **Set up Nextclade dataset:**

```bash
python scripts/setup_nextclade_dataset.py --email your.email@example.com
```

4. **Run the Nextclade pipeline:**

```bash
nextstrain build --docker .
```

5. **View the results:**

```bash
nextstrain view results/auspice/
```

# Automated Updates

To set up automated updates:

## 1. On Linux/macOS:

### Create a cron job:

```bash
# Edit crontab
crontab -e

# Add a line to run weekly (every Monday at 2 AM)
0 2 * * 1 cd /path/to/rvf-nextstrain && python scripts/download_ncbi_rvf_data.py --email your.email@example.com && nextstrain build . --docker
```

## 2. On Windows:

Use Windows Task Scheduler:

1. Create a new task
2. Set trigger (e.g., weekly)
3. Add action to run a batch script that executes:

```bash
cd C:\path\to\rvf-nextstrain
python scripts\download_ncbi_rvf_data.py --email your.email@example.com
nextstrain build .
```

# Project Structure

- config: Configuration files
- data: Raw sequences and metadata
- scripts: Python scripts for data processing
- workflow: Snakemake rules
- results: Pipeline outputs
- datasets: Nextclade reference data

# Customization

- To customize for different segments or parameters, modify `config.yaml`.
