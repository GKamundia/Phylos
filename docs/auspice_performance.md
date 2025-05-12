# RVF-Nextstrain Visualization Performance Optimization

The RVF-Nextstrain dashboard is optimized for performance through careful subsampling strategies. Large phylogenetic trees can cause browser performance issues in Auspice, especially when dealing with thousands of sequences.

## Current Subsampling Strategy

Our current approach to manage dataset size for optimal visualization performance:

1. **Group-based Subsampling**:
   - Sequences are grouped by `country year` to ensure geographic and temporal representation
   - Maximum sequences per build is configured in `config/master_config.yaml` as 300 sequences
   - Priority is given to sequences from human hosts and major livestock species

```yaml
# From pathogens/rvf/config/config.yaml
subsample:
  max_sequences: 300
  group_by: "country year"
  priorities:
    host:
      type: "categorical"
      focus:
        - "Homo sapiens"
        - "Ovis aries"
        - "Bos taurus"
```
