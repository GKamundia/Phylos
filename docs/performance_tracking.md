# Performance Tracking in RVF-Nextstrain

This document explains how performance tracking works in the RVF-Nextstrain pipeline and how to interpret the results.

## Built-in Snakemake Benchmarking

The pipeline automatically tracks performance metrics for each rule using Snakemake's benchmarking feature. When a rule runs, it generates a benchmark file in the `benchmarks/` directory with detailed metrics:

- Wall-clock runtime
- CPU time
- Memory usage (RSS, VSS, USS, PSS)
- I/O statistics
- CPU load

Example benchmark directive in a rule:

```smk
benchmark:
    f"benchmarks/tree_{output_prefix}.txt"
```
