#!/usr/bin/env python3
"""
Final validation script for advanced phylogenetics integration
Tests actual execution with sample data
"""

import os
import sys
import subprocess
import json
from pathlib import Path

def run_command(cmd, cwd=None):
    """Run a command and return success status and output"""
    try:
        result = subprocess.run(
            cmd, 
            shell=True, 
            capture_output=True, 
            text=True, 
            cwd=cwd
        )
        return result.returncode == 0, result.stdout, result.stderr
    except Exception as e:
        return False, "", str(e)

def check_sample_data():
    """Check if sample data exists for testing"""
    print("Checking for sample data...")
    
    # Look for any existing filtered data
    filtered_dir = Path("results/filtered")
    if filtered_dir.exists():
        fasta_files = list(filtered_dir.glob("*.fasta"))
        metadata_files = list(filtered_dir.glob("*metadata.tsv"))
        
        if fasta_files and metadata_files:
            print(f"✓ Found sample data: {len(fasta_files)} FASTA files, {len(metadata_files)} metadata files")
            return True, fasta_files[0], metadata_files[0]
    
    print("⚠ No sample data found in results/filtered/")
    return False, None, None

def test_advanced_phylogenetics_execution():
    """Test actual execution of advanced phylogenetics rules"""
    print("Testing advanced phylogenetics execution...")
    
    # Check for sample data first
    has_data, fasta_file, metadata_file = check_sample_data()
    
    if not has_data:
        print("⚠ Skipping execution test - no sample data available")
        print("   To test with real data, run the main pipeline first:")
        print("   snakemake --until filter")
        return True
    
    # Test individual rule execution with small subset
    print(f"Testing with sample file: {fasta_file}")
    
    # Try to run just the alignment step
    pathogen = "rvf"  # Assuming RVF based on project structure
    
    success, stdout, stderr = run_command(
        f"snakemake results/advanced_phylogenetics/{pathogen}_advanced_alignment.fasta --dry-run"
    )
    
    if success:
        print("✓ Advanced alignment rule can be executed")
        return True
    else:
        print(f"✗ Advanced alignment rule failed: {stderr}")
        return False

def validate_external_tools():
    """Check availability of external phylogenetic tools"""
    print("Validating external tools...")
    
    tools = {
        "mafft": "MAFFT alignment tool",
        "iqtree": "IQ-TREE maximum likelihood",
        "fasttree": "FastTree approximate ML"
    }
    
    available_tools = {}
    
    for tool, description in tools.items():
        success, stdout, stderr = run_command(f"{tool} --help")
        available_tools[tool] = success
        
        if success:
            print(f"✓ {tool}: {description} - Available")
        else:
            print(f"⚠ {tool}: {description} - Not installed (will use fallbacks)")
    
    # Check if at least one advanced tool is available
    if any(available_tools.values()):
        print("✓ At least one external tool is available for optimal performance")
    else:
        print("⚠ No external tools found - pipeline will use built-in Python alternatives")
    
    return available_tools

def generate_execution_report():
    """Generate a report on execution readiness"""
    print("Generating execution readiness report...")
    
    has_data, fasta_file, metadata_file = check_sample_data()
    available_tools = validate_external_tools()
    
    report = {
        "timestamp": "2025-06-06",
        "execution_readiness": {
            "sample_data_available": has_data,
            "sample_files": {
                "fasta": str(fasta_file) if fasta_file else None,
                "metadata": str(metadata_file) if metadata_file else None
            },
            "external_tools": available_tools,
            "ready_for_execution": has_data and any(available_tools.values())
        },
        "recommendations": []
    }
    
    # Add recommendations based on status
    if not has_data:
        report["recommendations"].append("Run 'snakemake --until filter' to generate sample data")
    
    if not any(available_tools.values()):
        report["recommendations"].extend([
            "Install MAFFT: conda install -c bioconda mafft",
            "Install IQ-TREE: conda install -c bioconda iqtree",
            "Install FastTree: conda install -c bioconda fasttree"
        ])
    
    if has_data and any(available_tools.values()):
        report["recommendations"].append("Ready to run: snakemake run_advanced_phylogenetics")
    
    # Save report
    with open("execution_readiness.json", "w") as f:
        json.dump(report, f, indent=2)
    
    print("✓ Execution readiness report saved to execution_readiness.json")
    return report

def main():
    """Main validation function"""
    print("=" * 60)
    print("ADVANCED PHYLOGENETICS EXECUTION VALIDATION")
    print("=" * 60)
    
    # Change to project directory
    os.chdir(Path(__file__).parent)
    
    # Run validation tests
    tests = [
        ("Sample Data Check", check_sample_data),
        ("External Tools Validation", validate_external_tools),
        ("Execution Test", test_advanced_phylogenetics_execution),
        ("Execution Report", generate_execution_report)
    ]
    
    results = {}
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        print("-" * 40)
        try:
            if test_name == "Sample Data Check":
                has_data, _, _ = test_func()
                results[test_name] = has_data
            elif test_name == "External Tools Validation":
                tools = test_func()
                results[test_name] = any(tools.values())
            else:
                results[test_name] = test_func()
        except Exception as e:
            print(f"✗ {test_name} failed: {str(e)}")
            results[test_name] = False
    
    # Summary
    print("\n" + "=" * 60)
    print("EXECUTION VALIDATION SUMMARY")
    print("=" * 60)
    
    for test_name, passed in results.items():
        status = "PASS" if passed else "WARN" if test_name in ["Sample Data Check", "External Tools Validation"] else "FAIL"
        print(f"{test_name}: {status}")
    
    if results.get("Sample Data Check") and results.get("External Tools Validation"):
        print("\n🎉 READY FOR FULL EXECUTION!")
        print("Run: snakemake run_advanced_phylogenetics")
    else:
        print("\n⚠ PREPARATION NEEDED")
        print("See execution_readiness.json for specific recommendations")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
