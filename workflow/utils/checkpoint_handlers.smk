"""
Handlers for Snakemake checkpoints
"""
import json

def get_available_segments(wildcards):
    """
    Dynamic function to get available segments after the split_by_segment checkpoint
    
    This allows the workflow to adapt if some segments don't have enough data
    """
    # Get the split_by_segment checkpoint output
    checkpoint_output = checkpoints.split_by_segment.get(**wildcards).output
    
    # Find which segments were successfully processed
    available_segments = []
    for segment in segments:
        segment_file = f"results/segments/{segment}/raw/{wildcards.prefix}_{segment}_sequences.fasta"
        if os.path.exists(segment_file) and os.path.getsize(segment_file) > 0:
            available_segments.append(segment)
    
    if not available_segments:
        print(f"Warning: No valid segment data found for {wildcards.prefix}")
        # Default to attempting all segments (will fail gracefully later)
        return segments
    
    print(f"Found data for segments: {', '.join(available_segments)}")
    return available_segments

def get_metadata_result(wildcards):
    """
    Dynamic function to handle conditional metadata validation checkpoint
    """
    metadata_output = checkpoints.prepare_metadata.get(**wildcards).output
    
    # Check validation report for strict mode decision
    if config.get("workflow", {}).get("strict_metadata", False):
        validation_file = metadata_output.report  # Fixed: renamed from report_file
        try:
            with open(validation_file) as f:
                report_data = json.load(f)  # Renamed from 'report' to 'report_data'
                if report_data.get("invalid_records", 0) > 0:
                    print(f"Warning: {report_data.get('invalid_records')} invalid metadata records found")
                    if config.get("workflow", {}).get("fail_on_invalid", False):
                        raise ValueError(f"Strict metadata validation failed. See {validation_file} for details")
        except Exception as e:
            print(f"Warning: Could not validate metadata report: {e}")
    
    return metadata_output.metadata