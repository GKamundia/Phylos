#!/usr/bin/env python3
"""
Set up a basic Nextclade dataset for Rift Valley fever virus
"""

import os
import argparse
import sys
from Bio import Entrez, SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Set up Nextclade dataset for RVF")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--output-dir", default="nextclade/datasets/rvf",
                        help="Output directory for Nextclade dataset")
    parser.add_argument("--reference-accession", default="DQ375404",
                        help="Accession for L segment reference")
    return parser.parse_args()

def download_reference(accession, email):
    """Download reference genome from NCBI"""
    Entrez.email = email
    
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, 
                              rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error downloading reference: {e}")
        return None

def create_genemap(record):
    """Create a simple GFF from GenBank features"""
    gff_lines = ["##gff-version 3"]
    
    seq_id = record.id
    for feature in record.features:
        if feature.type == "CDS" or feature.type == "gene":
            # Get basic feature info
            start = feature.location.start + 1  # GFF is 1-based
            end = feature.location.end
            strand = '+' if feature.location.strand == 1 else '-'
            
            # Get the feature name/ID
            feature_id = ""
            if "gene" in feature.qualifiers:
                feature_id = feature.qualifiers["gene"][0]
            elif "product" in feature.qualifiers:
                feature_id = feature.qualifiers["product"][0]
            else:
                feature_id = f"{feature.type}_{start}_{end}"
            
            # Format GFF line
            attributes = f"ID={feature_id};Name={feature_id}"
            if "product" in feature.qualifiers:
                attributes += f";product={feature.qualifiers['product'][0]}"
            
            gff_line = f"{seq_id}\tGenBank\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}"
            gff_lines.append(gff_line)
    
    return "\n".join(gff_lines)

def create_qc_config():
    """Create a basic QC config for RVF"""
    qc_json = """{
  "schemaVersion": "1.3.0",
  "privateMutations": {
    "enabled": true,
    "typical": 0,
    "cutoff": 15,
    "weightEmpty": 6000
  },
  "missingData": {
    "enabled": true,
    "typical": 0,
    "cutoff": 2700,
    "weightEmpty": 6000
  },
  "snpClusters": {
    "enabled": true,
    "windowSize": 100,
    "clusterCutOff": 6,
    "scoreWeight": 50
  },
  "frameShifts": {
    "enabled": false
  },
  "stopCodons": {
    "enabled": false
  }
}"""
    return qc_json

def main():
    args = parse_args()
    
    # Make sure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Download reference genome
    print(f"Downloading reference genome (accession: {args.reference_accession})...")
    reference = download_reference(args.reference_accession, args.email)
    
    if not reference:
        print("Failed to download reference genome")
        return 1
    
    # Write reference FASTA
    ref_path = os.path.join(args.output_dir, "reference.fasta")
    with open(ref_path, 'w') as f:
        SeqIO.write(reference, f, "fasta")
    print(f"Wrote reference genome to {ref_path}")
    
    # Create and write genemap.gff
    genemap = create_genemap(reference)
    genemap_path = os.path.join(args.output_dir, "genemap.gff")
    with open(genemap_path, 'w') as f:
        f.write(genemap)
    print(f"Created gene map at {genemap_path}")
    
    # Create and write QC config
    qc_config = create_qc_config()
    qc_path = os.path.join(args.output_dir, "qc.json")
    with open(qc_path, 'w') as f:
        f.write(qc_config)
    print(f"Created QC config at {qc_path}")
    
    print(f"Nextclade dataset for RVF set up in {args.output_dir}")
    print("Next steps: Review and adjust files if needed")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())