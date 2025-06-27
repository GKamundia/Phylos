# Filter System Update - Removal of Generic min_length Configuration

## Summary of Changes

The filtering system has been updated to remove generic `min_length` configuration-based filtering while retaining the exact reference length filtering for RVF segments.

## Changes Made

### 1. Updated `scripts/run_filter.py`

- **Added REFERENCE_LENGTHS definition**: Hardcoded the exact reference lengths for RVF segments
  ```python
  REFERENCE_LENGTHS = {
      'L': 6404,  # L segment exact reference length
      'M': 3885,  # M segment exact reference length
      'S': 1690   # S segment exact reference length
  }
  ```
- **Updated documentation**: Clarified that filtering uses exact reference lengths, not generic min_length
- **Kept exact reference length filtering**: The script still filters sequences to match exact reference lengths

### 2. Updated `pathogens/templates/config/config.yaml`

- **Removed min_length parameter**: No longer part of the template configuration
- **Added documentation note**: Explains that exact reference length filtering is handled in the script

## Current Filtering Behavior

The filtering system now works as follows:

1. **Segment identification**: Based on 'Segment' column in metadata or sequence naming patterns
2. **Completeness filtering**: Only sequences with `Nuc_Completeness == 'complete'`
3. **Exact reference length filtering**: Sequences must match exact reference lengths:
   - L segment: exactly 6404 bp
   - M segment: exactly 3885 bp
   - S segment: exactly 1690 bp
4. **Metadata filtering**: Sequences must have non-empty country and date fields

## What Was Removed

- Generic `min_length` configuration parameter from pathogen configs
- Flexible minimum length thresholds that could be set in configuration
- Any filtering rules that used configurable minimum lengths

## What Was Kept

- **Exact reference length filtering**: This ensures phylogenetic analysis quality
- **All other filtering criteria**: Completeness, country, date validation
- **QC min_length (200bp)**: General QC validation threshold remains unchanged

## Benefits

1. **Simplified configuration**: No confusing min_length parameters in configs
2. **Consistent quality**: All sequences match exact reference lengths
3. **Clear filtering logic**: Segment-specific exact matching vs generic thresholds
4. **Better phylogenetic analysis**: Uniform sequence lengths improve alignment quality

## Files Modified

- `scripts/run_filter.py` - Added REFERENCE_LENGTHS definition and updated docs
- `pathogens/templates/config/config.yaml` - Removed min_length parameter

## Files NOT Modified

- `config/qc_config.json` - QC min_length (200bp) kept for general validation
- `scripts/split_by_segment.py` - Already uses exact reference length filtering
- Workflow rules - No changes needed since min_length params weren't being used

## Testing

- Both modified scripts compile without syntax errors
- Exact reference length filtering logic remains intact
- No breaking changes to the filtering pipeline
