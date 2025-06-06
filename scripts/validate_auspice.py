#!/usr/bin/env python3
"""
Validate the generated Auspice JSON file
"""

import json
import os

def validate_auspice_json():
    json_file = "auspice/rvf.json"
    
    if not os.path.exists(json_file):
        print(f"Error: {json_file} not found")
        return False
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        print("🎉 Auspice JSON Validation Results:")
        print("="*50)
        
        # Check version
        print(f"Version: {data.get('version', 'Unknown')}")
        
        # Check meta information
        meta = data.get('meta', {})
        print(f"Title: {meta.get('title', 'Unknown')}")
        print(f"Updated: {meta.get('updated', 'Unknown')}")
        
        # Check colorings
        colorings = meta.get('colorings', [])
        print(f"Available colorings: {len(colorings)}")
        for coloring in colorings:
            print(f"  - {coloring.get('key', 'Unknown')}: {coloring.get('title', 'Unknown')}")
        
        # Check tree structure
        tree = data.get('tree', {})
        children = tree.get('children', [])
        print(f"Tree children: {len(children)}")
        
        # Check nodes
        nodes = data.get('nodes', {})
        print(f"Total nodes: {len(nodes)}")
        
        # Check sample node attributes
        if nodes:
            sample_node = list(nodes.values())[0]
            node_attrs = sample_node.get('node_attrs', {})
            print(f"Sample node attributes: {list(node_attrs.keys())}")
        
        # Check geo resolutions
        geo_resolutions = meta.get('geo_resolutions', [])
        if geo_resolutions:
            demes = geo_resolutions[0].get('demes', {})
            print(f"Countries in geo resolution: {len(demes)}")
            if demes:
                print(f"Sample countries: {list(demes.keys())[:5]}")
        
        # File size
        file_size = os.path.getsize(json_file)
        print(f"File size: {file_size:,} bytes ({file_size/1024/1024:.1f} MB)")
        
        print("="*50)
        print("✅ JSON validation successful!")
        return True
        
    except Exception as e:
        print(f"❌ Error validating JSON: {e}")
        return False

if __name__ == "__main__":
    validate_auspice_json()
