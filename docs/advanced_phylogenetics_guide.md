# Advanced Phylogenetic Analysis Guide

This guide covers the advanced phylogenetic analysis capabilities implemented for the RVF Nextstrain project. These tools provide enhanced sequence alignment, tree building, and temporal analysis features beyond the standard Nextstrain pipeline.

## Overview

The advanced phylogenetic analysis system consists of three main components:

1. **Advanced Alignment** (`scripts/advanced_alignment.py`) - Multi-algorithm sequence alignment with fallback mechanisms
2. **Advanced Phylogenetics** (`scripts/advanced_phylogenetics.py`) - Multiple tree building methods with temporal integration
3. **Temporal Analysis** (`scripts/temporal_analysis.py`) - Comprehensive temporal pattern analysis and visualization

## Features

### Advanced Alignment

- **Multiple Algorithms**: MAFFT, MUSCLE, ClustalW, with automatic method selection
- **Fallback Mechanisms**: Built-in BioPython alignment when external tools are unavailable
- **Segmented Genome Support**: Specialized alignment for segmented viruses (L, M, S segments)
- **Quality Statistics**: Comprehensive alignment quality metrics
- **Reference-Guided Alignment**: Optional reference sequence support

### Advanced Phylogenetics

- **Multiple Tree Methods**: IQ-TREE, FastTree, RAxML, Neighbor Joining, Maximum Parsimony
- **Automatic Method Selection**: Chooses optimal method based on dataset characteristics
- **Temporal Integration**: Built-in temporal analysis capabilities
- **Comprehensive Statistics**: Tree depth, branch lengths, node counts
- **Model Selection**: Automatic or manual evolutionary model selection

### Temporal Analysis

- **Comprehensive Pattern Analysis**: Yearly, monthly, seasonal distributions
- **Geographic Temporal Patterns**: Analysis by country and division
- **Segment-Specific Analysis**: Temporal patterns for segmented viruses
- **Phylogenetic Temporal Signal**: Correlation between genetic and temporal distance
- **Rich Visualizations**: Multiple plot types for temporal patterns
- **Automated Reporting**: JSON results and human-readable summaries

## Installation

### Python Dependencies

The required Python packages are listed in `requirements.txt`:

```bash
pip install -r requirements.txt
```

### External Tools (Optional)

For enhanced performance, install these external tools:

- **MAFFT**: Fast multiple sequence alignment
- **IQ-TREE**: Maximum likelihood phylogenetic inference
- **FastTree**: Approximate maximum likelihood trees
- **RAxML**: Maximum likelihood phylogenetic inference
- **TreeTime**: Molecular clock analysis

## Usage Examples

### 1. Basic Advanced Alignment

```bash
# Auto-select alignment method
python scripts/advanced_alignment.py sequences.fasta aligned.fasta --method auto

# Use specific method with statistics
python scripts/advanced_alignment.py sequences.fasta aligned.fasta \
    --method mafft \
    --threads 8 \
    --stats-output alignment_stats.json
```

### 2. Segmented Virus Alignment

```bash
# Align segments separately
python scripts/advanced_alignment.py sequences.fasta output_dir/ \
    --metadata metadata.tsv \
    --segment-mode \
    --segments L M S \
    --stats-output segment_stats.json
```

### 3. Advanced Tree Building

```bash
# Auto-select tree method
python scripts/advanced_phylogenetics.py aligned.fasta tree.nwk \
    --method auto \
    --metadata metadata.tsv \
    --stats-output tree_stats.json

# Use specific method with temporal analysis
python scripts/advanced_phylogenetics.py aligned.fasta tree.nwk \
    --method iqtree \
    --model GTR+G \
    --temporal-analysis \
    --temporal-output-dir temporal_results/
```

### 4. Comprehensive Temporal Analysis

```bash
# Basic temporal analysis
python scripts/temporal_analysis.py metadata.tsv \
    --output-dir temporal_analysis/

# With phylogenetic tree
python scripts/temporal_analysis.py metadata.tsv \
    --tree tree.nwk \
    --output-dir temporal_analysis/

# Skip plot generation (faster)
python scripts/temporal_analysis.py metadata.tsv \
    --output-dir temporal_analysis/ \
    --no-plots
```

## Snakemake Integration

The advanced features are integrated into the Snakemake workflow via `workflow/core_rules/advanced_phylogenetics.smk`:

### Key Rules:

- `advanced_align`: Enhanced alignment with multiple algorithms
- `segment_align`: Segment-specific alignment for segmented viruses
- `advanced_tree`: Tree building with automatic method selection
- `temporal_analysis`: Comprehensive temporal pattern analysis
- `molecular_clock`: TreeTime-based molecular clock estimation
- `segment_phylogenetics`: Segment-specific phylogenetic analysis
- `phylogenetic_qc`: Quality control assessment
- `advanced_phylogenetic_summary`: Results aggregation

### Configuration

Add these options to your Snakemake config:

```yaml
# Advanced phylogenetic settings
align:
  method: "auto" # or "mafft", "muscle", "clustalw"

tree:
  method: "auto" # or "iqtree", "fasttree", "raxml", "neighbor_joining", "maximum_parsimony"
  model: "auto" # or specific model like "GTR+G"

temporal:
  clock_filter: 3

segments: ["L", "M", "S"] # For segmented viruses

resources:
  align:
    threads: 4
    mem_mb: 4000
  tree:
    threads: 8
    mem_mb: 8000
    runtime: 120
  temporal:
    threads: 2
    mem_mb: 4000
```

## Input Requirements

### Metadata Format

The temporal analysis expects metadata with date information. Supported date columns:

- `date`
- `Date`
- `Release_Date`
- `collection_date`
- `Collection_Date`

Supported date formats:

- YYYY-MM-DD (ISO format)
- DD-MMM-YYYY (e.g., "23-MAR-2025")
- DD/MM/YYYY or MM/DD/YYYY
- YYYY-MM or YYYY

Additional useful columns:

- `country`: Geographic analysis
- `division`: Sub-regional analysis
- `segment`: Segment-specific analysis for segmented viruses
- `strain`: Strain identifier for phylogenetic analysis

### Sequence Requirements

- FASTA format sequences
- Minimum 2 sequences for alignment
- Sequences should be of the same genomic segment/gene
- For best results, sequences should be pre-filtered for quality

## Output Files

### Alignment Output

- `*_aligned.fasta`: Aligned sequences
- `*_alignment_stats.json`: Alignment quality statistics

### Phylogenetic Output

- `*.nwk`: Newick format phylogenetic tree
- `*_tree_stats.json`: Tree statistics and metadata
- `temporal_analysis/`: Temporal analysis results (if metadata provided)

### Temporal Analysis Output

- `temporal_analysis_results.json`: Complete analysis results
- `temporal_analysis_summary.txt`: Human-readable summary
- `yearly_distribution.png`: Yearly sampling plots
- `monthly_patterns.png`: Monthly and seasonal patterns
- `geographic_temporal.png`: Geographic temporal patterns
- `segment_temporal.png`: Segment-specific patterns (if applicable)

## Performance Considerations

### Method Selection Guidelines:

- **Small datasets (<50 sequences)**: IQ-TREE for best accuracy
- **Medium datasets (50-200 sequences)**: FastTree for speed-accuracy balance
- **Large datasets (>200 sequences)**: Neighbor Joining for speed

### Resource Requirements:

- **Alignment**: 4-8 GB RAM for typical datasets
- **Tree Building**: 8-16 GB RAM, scales with sequence number and length
- **Temporal Analysis**: 2-4 GB RAM, mainly I/O bound

### Optimization Tips:

1. Use `--threads` parameter to utilize multiple CPU cores
2. For large datasets, consider subsampling before analysis
3. Use `--no-plots` for temporal analysis when visualizations aren't needed
4. External tools (MAFFT, IQ-TREE) are faster than BioPython fallbacks

## Troubleshooting

### Common Issues:

1. **External tools not found**: Scripts will automatically fall back to BioPython implementations
2. **Date parsing errors**: Check date format and column names in metadata
3. **Memory errors**: Reduce dataset size or increase available memory
4. **Empty output**: Check input file formats and paths

### Debug Mode:

Add `--verbose` or check log files in the `logs/` directory for detailed error information.

## Integration with Main Workflow

To integrate these features into your main Snakefile:

```python
# Include the advanced phylogenetic rules
include: "workflow/core_rules/advanced_phylogenetics.smk"

# Add advanced targets to your main rule
rule all:
    input:
        # ... existing targets ...
        f"results/advanced_phylogenetic_summary_{output_prefix}.json",
        f"results/temporal/{output_prefix}/temporal_analysis_summary.txt"
```

## Future Enhancements

Planned improvements include:

- Integration with more external phylogenetic tools
- Enhanced molecular clock models
- Automated outbreak detection
- Real-time temporal monitoring
- Integration with geographic mapping tools

## Support

For issues or questions:

1. Check the log files in `logs/` directory
2. Verify input file formats match requirements
3. Ensure all dependencies are properly installed
4. Review this documentation for configuration options

The advanced phylogenetic analysis system is designed to be robust and provide fallback options when external tools are unavailable, ensuring reliable operation across different computing environments.
