FROM nextstrain/base:build-20250505T154520Z

# Set working directory
WORKDIR /nextstrain/rvf-nextstrain

# Install additional Python dependencies with pinned versions
RUN pip install --no-cache-dir \
    biopython==1.82 \
    pandas==2.1.2 \
    pyyaml==6.0.1 \
    jsonschema==4.19.1 \
    pytest==7.4.3 \
    requests==2.31.0

# Set environment variables 
ENV PYTHONPATH="/nextstrain/rvf-nextstrain"
ENV PATH="/nextstrain/rvf-nextstrain/scripts:${PATH:-}"

# Set default command
CMD ["/bin/bash"]