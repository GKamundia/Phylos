#!/usr/bin/env python3
"""
Test script to validate advanced phylogenetics integration
"""

import os
import sys
import yaml
import json
from pathlib import Path

# Add scripts directory to path
sys.path.append('scripts')

def test_configuration_loading():
    """Test that the advanced phylogenetics configuration loads correctly"""
    print("Testing configuration loading...")
    
    try:
        # Load master config
        with open('config/master_config.yaml', 'r') as f:
            master_config = yaml.safe_load(f)
        
        # Load pathogen-specific config
        with open('pathogens/rvf/config/config.yaml', 'r') as f:
            rvf_config = yaml.safe_load(f)
        
        # Check if advanced phylogenetics is enabled
        master_advanced = master_config.get('common', {}).get('advanced_phylogenetics', {})
        rvf_advanced = rvf_config.get('advanced_phylogenetics', {})
        
        print(f"Master config has advanced phylogenetics: {bool(master_advanced)}")
        print(f"RVF config has advanced phylogenetics: {bool(rvf_advanced)}")
        print(f"Advanced phylogenetics enabled: {rvf_advanced.get('enabled', False)}")
        
        # Merge configurations as done in Snakefile
        merged_config = {**master_advanced, **rvf_advanced}
        print(f"Merged configuration: {merged_config.get('enabled', False)}")
        
        return True
        
    except Exception as e:
        print(f"Configuration loading failed: {e}")
        return False

def test_advanced_modules():
    """Test that advanced phylogenetics modules can be imported"""
    print("\nTesting module imports...")
    
    try:
        from advanced_alignment import AdvancedAligner
        from advanced_phylogenetics import AdvancedPhylogeneticAnalyzer
        from temporal_analysis import TemporalAnalysisEngine
        
        print("✓ All advanced modules imported successfully")
        return True
        
    except ImportError as e:
        print(f"✗ Module import failed: {e}")
        return False

def test_directory_structure():
    """Test that required directories exist"""
    print("\nTesting directory structure...")
    
    required_dirs = [
        "results/advanced_phylogenetics",
        "workflow/core_rules",
        "scripts"
    ]
    
    all_exist = True
    for directory in required_dirs:
        if os.path.exists(directory):
            print(f"✓ {directory} exists")
        else:
            print(f"✗ {directory} missing")
            all_exist = False
    
    # Check for advanced phylogenetics rule file
    rule_file = "workflow/core_rules/advanced_phylogenetics.smk"
    if os.path.exists(rule_file):
        print(f"✓ {rule_file} exists")
    else:
        print(f"✗ {rule_file} missing")
        all_exist = False
    
    return all_exist

def test_data_availability():
    """Test that test data is available"""
    print("\nTesting data availability...")
    
    test_files = [
        "data/metadata/rvf_metadata_with_segments.tsv",
        "data/sequences/raw/rvf_sequences.fasta"
    ]
    
    available = True
    for file_path in test_files:
        if os.path.exists(file_path):
            print(f"✓ {file_path} available")
        else:
            print(f"✗ {file_path} not found")
            available = False
    
    return available

def test_basic_functionality():
    """Test basic functionality of advanced phylogenetics modules"""
    print("\nTesting basic functionality...")
    
    try:
        # Test AdvancedAligner
        from advanced_alignment import AdvancedAligner
        aligner = AdvancedAligner()
        print("✓ AdvancedAligner initialized")
        
        # Test AdvancedPhylogeneticAnalyzer
        from advanced_phylogenetics import AdvancedPhylogeneticAnalyzer
        analyzer = AdvancedPhylogeneticAnalyzer()
        print("✓ AdvancedPhylogeneticAnalyzer initialized")
        
        # Test TemporalAnalysisEngine
        from temporal_analysis import TemporalAnalysisEngine
        temporal = TemporalAnalysisEngine()
        print("✓ TemporalAnalysisEngine initialized")
        
        return True
        
    except Exception as e:
        print(f"✗ Functionality test failed: {e}")
        return False

def main():
    """Run all integration tests"""
    print("=" * 60)
    print("ADVANCED PHYLOGENETICS INTEGRATION TEST")
    print("=" * 60)
    
    tests = [
        ("Configuration Loading", test_configuration_loading),
        ("Module Imports", test_advanced_modules),
        ("Directory Structure", test_directory_structure),
        ("Data Availability", test_data_availability),
        ("Basic Functionality", test_basic_functionality)
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        print("-" * 40)
        result = test_func()
        results.append((test_name, result))
    
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    all_passed = True
    for test_name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"{test_name}: {status}")
        if not result:
            all_passed = False
    
    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 ALL TESTS PASSED - Integration successful!")
        print("Advanced phylogenetics capabilities are ready to use.")
    else:
        print("❌ SOME TESTS FAILED - Integration needs attention.")
        print("Please check the failed components before proceeding.")
    
    print("=" * 60)
    
    return all_passed

if __name__ == "__main__":
    # Change to the correct directory
    os.chdir(Path(__file__).parent)
    success = main()
    sys.exit(0 if success else 1)
