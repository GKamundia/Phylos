# Base image with specific version for reproducibility
FROM nextstrain/base:20.11.1

LABEL maintainer="Kamundia <gkamush50@gmail.com>"
LABEL description="Reproducible Nextstrain environment for RVF analysis"
LABEL version="1.0.0"

# Set working directory
WORKDIR /nextstrain/rvf-nextstrain

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    git \
    vim \
    parallel \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies with pinned versions
RUN pip install --no-cache-dir \
    biopython==1.82 \
    pandas==2.1.2 \
    pyarrow==15.0.0 \
    pyyaml==6.0.1 \
    jsonschema==4.19.1 \
    pytest==7.4.3 \
    requests==2.31.0 \
    matplotlib==3.8.0 \
    seaborn==0.13.0 \
    geopy==2.4.1 \
    nextclade==2.14.0

# Install latest Nextstrain CLI for backwards compatibility
RUN pip install nextstrain-cli

# Create necessary directories
RUN mkdir -p /nextstrain/rvf-nextstrain/data \
    /nextstrain/rvf-nextstrain/config \
    /nextstrain/rvf-nextstrain/results \
    /nextstrain/rvf-nextstrain/logs \
    /nextstrain/rvf-nextstrain/benchmarks

# Set environment variables
ENV PYTHONPATH="/nextstrain/rvf-nextstrain"
ENV PYTHONPATH="${PYTHONPATH}:/nextstrain/rvf-nextstrain"
ENV PATH="/nextstrain/rvf-nextstrain/scripts:${PATH}"

# Set up entry point for Nextstrain CLI
ENTRYPOINT ["nextstrain"]
CMD ["--help"]