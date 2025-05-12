#!/usr/bin/env python3
"""
Standalone metadata validation tool for RVF Nextstrain analysis
"""

import os
import sys
import json
import pandas as pd
import argparse
import logging
import jsonschema
from jsonschema import validate, ValidationError

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate metadata against JSON schema")
    parser.add_argument('metadata', help="Path to metadata TSV file")
    parser.add_argument('--schema', default="config/metadata_schema.json", 
                        help="Path to JSON schema file")
    parser.add_argument('--report', default=None,
                        help="Path to write validation report (optional)")
    return parser.parse_args()

def load_schema(schema_file):
    """Load and parse the JSON schema file"""
    try:
        with open(schema_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load schema file: {e}")
        sys.exit(1)

def validate_record(record, schema):
    """Validate a single metadata record against the schema"""
    try:
        validate(instance=record, schema=schema)
        return True, None
    except ValidationError as e:
        return False, str(e)

def main():
    args = parse_args()
    
    # Load schema
    schema = load_schema(args.schema)
    logger.info(f"Loaded schema from {args.schema}")
    
    # Load metadata
    logger.info(f"Loading metadata from {args.metadata}")
    try:
        metadata = pd.read_csv(args.metadata, sep='\t', dtype=str)
    except Exception as e:
        logger.error(f"Error loading metadata: {e}")
        return 1
    
    logger.info(f"Loaded {len(metadata)} records")
    
    # Validate each record
    valid_count = 0
    invalid_records = []
    
    for idx, row in metadata.iterrows():
        record = row.to_dict()
        valid, error = validate_record(record, schema)
        
        if valid:
            valid_count += 1
        else:
            invalid_records.append({
                'row': int(idx),
                'strain': row.get('strain', 'UNKNOWN'),
                'error': error
            })
    
    # Generate summary
    logger.info(f"Valid records: {valid_count}/{len(metadata)} ({valid_count/len(metadata)*100:.1f}%)")
    if invalid_records:
        logger.warning(f"Invalid records: {len(invalid_records)}")
        
        # Write detailed report if requested
        if args.report:
            report = {
                'summary': {
                    'total_records': len(metadata),
                    'valid_records': valid_count,
                    'invalid_records': len(invalid_records),
                    'validation_rate': f"{valid_count/len(metadata)*100:.1f}%"
                },
                'invalid_records': invalid_records
            }
            
            with open(args.report, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Validation report written to {args.report}")
    else:
        logger.info("All records passed validation")
    
    return 0 if len(invalid_records) == 0 else 1

if __name__ == "__main__":
    sys.exit(main())