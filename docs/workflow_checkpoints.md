# Workflow Checkpoints

This document describes the critical checkpoints in the RVF-Nextstrain pipeline and how they're used to make the workflow more resilient.

## Implemented Checkpoints

### 1. Metadata Preparation (`prepare_metadata`)

**Purpose**: Validates and standardizes metadata before proceeding with analysis.

**Benefits**:

- Allows the workflow to fail early if metadata doesn't meet requirements
- In strict mode, prevents invalid data from entering the pipeline
- Creates a validation report for quality assurance

**Downstream Effects**:

- The `get_metadata_result()` function uses this checkpoint to determine if analysis should continue
- Configurable via `workflow.strict_metadata` and `workflow.fail_on_invalid` parameters

### 2. Segment Splitting (`split_by_segment`)

**Purpose**: In multi-segment mode, determines which genome segments have sufficient data.

**Benefits**:

- Makes the workflow adaptable to available data
- Prevents failures when some segments have insufficient sequences
- Enables partial analyses when not all segments are available

**Downstream Effects**:

- The `get_available_segments()` function uses this checkpoint to determine which segments to process
- Affects the segment-specific analysis paths and final output files

## Resuming from Checkpoints

To resume a failed workflow from the last successful checkpoint:

```bash
# Resume the entire workflow
nextstrain build . --restart-times 1

# Resume from a specific checkpoint
nextstrain build . --forcerun prepare_metadata
```
