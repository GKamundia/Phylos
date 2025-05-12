#!/usr/bin/env python3
"""
Cache manager for Auspice JSON files to improve dashboard performance
"""

import os
import json
import time
import hashlib
import logging
import shutil
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Union

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("cache_manager")

class AuspiceJsonCache:
    """
    Manages cached versions of Auspice JSON files
    """
    
    def __init__(self, 
                cache_dir: str = "cache/auspice",
                source_dir: str = "results/auspice",
                max_age_days: int = 30,
                max_versions: int = 10):
        """
        Initialize the cache manager
        
        Args:
            cache_dir: Directory to store cached files
            source_dir: Directory containing source Auspice JSON files
            max_age_days: Maximum age of cached files in days
            max_versions: Maximum number of versions to keep per file
        """
        self.cache_dir = Path(cache_dir)
        self.source_dir = Path(source_dir)
        self.max_age_days = max_age_days
        self.max_versions = max_versions
        
        # Create cache directory if it doesn't exist
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Cache index tracks available files and their versions
        self.cache_index_path = self.cache_dir / "cache_index.json"
        self.cache_index = self._load_cache_index()
    
    def _load_cache_index(self) -> Dict:
        """Load the cache index or create if it doesn't exist"""
        if self.cache_index_path.exists():
            try:
                with open(self.cache_index_path, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                logger.error(f"Failed to load cache index: {e}")
        
        # Default index structure
        return {
            "last_updated": datetime.now().isoformat(),
            "files": {}
        }
    
    def _save_cache_index(self):
        """Save the cache index"""
        self.cache_index["last_updated"] = datetime.now().isoformat()
        try:
            with open(self.cache_index_path, 'w') as f:
                json.dump(self.cache_index, f, indent=2)
        except IOError as e:
            logger.error(f"Failed to save cache index: {e}")
    
    def update_cache(self, force: bool = False):
        """
        Update cache with latest Auspice JSON files
        
        Args:
            force: If True, force update all files regardless of changes
        """
        # Find all JSON files in source directory
        source_files = list(self.source_dir.glob("*.json"))
        
        for source_file in source_files:
            file_name = source_file.name
            
            # Calculate file hash to detect changes
            file_hash = self._calculate_file_hash(source_file)
            
            # Check if file exists in cache and has changed
            file_entry = self.cache_index["files"].get(file_name, {"versions": []})
            latest_version = file_entry["versions"][0] if file_entry["versions"] else None
            
            if force or not latest_version or latest_version["hash"] != file_hash:
                # File is new or has changed, add to cache
                timestamp = datetime.now().isoformat()
                version_id = f"{int(time.time())}_{file_hash[:8]}"
                
                # Create cache file path with version
                cache_path = self.cache_dir / f"{file_name.replace('.json', '')}_{version_id}.json"
                
                # Copy file to cache
                try:
                    shutil.copy2(source_file, cache_path)
                    
                    # Create optimized versions
                    self._create_optimized_versions(source_file, cache_path, version_id)
                    
                    # Update index
                    new_version = {
                        "version_id": version_id,
                        "timestamp": timestamp,
                        "hash": file_hash,
                        "path": str(cache_path),
                        "size": os.path.getsize(source_file),
                        "optimized": {
                            "minimal": str(cache_path).replace(".json", "_minimal.json"),
                            "tree_only": str(cache_path).replace(".json", "_tree_only.json"),
                            "metadata_only": str(cache_path).replace(".json", "_metadata_only.json")
                        }
                    }
                    
                    file_entry["versions"].insert(0, new_version)
                    self.cache_index["files"][file_name] = file_entry
                    
                    logger.info(f"Added new cache version for {file_name}: {version_id}")
                except IOError as e:
                    logger.error(f"Failed to cache {file_name}: {e}")
        
        # Clean up old versions
        self._clean_cache()
        
        # Save updated index
        self._save_cache_index()
    
    def _create_optimized_versions(self, source_file: Path, cache_path: Path, version_id: str):
        """
        Create optimized versions of the cached file
        
        Args:
            source_file: Path to source file
            cache_path: Path to cached file
            version_id: Version identifier
        """
        try:
            # Load the full Auspice JSON
            with open(source_file, 'r') as f:
                data = json.load(f)
            
            # Create minimal version (no node data)
            minimal_path = str(cache_path).replace(".json", "_minimal.json")
            minimal_data = {k: v for k, v in data.items() if k != "tree"}
            if "tree" in data:
                # Include tree structure but remove most node attributes
                minimal_data["tree"] = self._simplify_tree(data["tree"])
            
            # Create tree-only version
            tree_path = str(cache_path).replace(".json", "_tree_only.json")
            tree_data = {"tree": data.get("tree", {}), "version": data.get("version", {})}
            
            # Create metadata-only version
            meta_path = str(cache_path).replace(".json", "_metadata_only.json")
            meta_data = {k: v for k, v in data.items() if k != "tree"}
            
            # Write optimized versions
            with open(minimal_path, 'w') as f:
                json.dump(minimal_data, f)
            
            with open(tree_path, 'w') as f:
                json.dump(tree_data, f)
                
            with open(meta_path, 'w') as f:
                json.dump(meta_data, f)
            
            logger.info(f"Created optimized versions for {source_file.name}")
        except Exception as e:
            logger.error(f"Failed to create optimized versions for {source_file.name}: {e}")
    
    def _simplify_tree(self, tree: Dict) -> Dict:
        """
        Simplify tree structure for minimal version
        
        Args:
            tree: Tree data structure
        
        Returns:
            Simplified tree
        """
        # Function to recursively simplify nodes
        def simplify_node(node):
            # Keep only essential node attributes
            essential_attrs = ["name", "children", "node_attrs"]
            simplified = {k: v for k, v in node.items() if k in essential_attrs}
            
            # Further simplify node_attrs to only include essential information
            if "node_attrs" in simplified:
                essential_node_attrs = ["div", "num_date", "country", "region"]
                simplified["node_attrs"] = {k: v for k, v in simplified["node_attrs"].items() 
                                          if k in essential_node_attrs}
            
            # Process children recursively
            if "children" in simplified:
                simplified["children"] = [simplify_node(child) for child in simplified["children"]]
            
            return simplified
        
        # Start with the root node
        return simplify_node(tree)
    
    def _calculate_file_hash(self, file_path: Path) -> str:
        """Calculate SHA-256 hash of a file"""
        hasher = hashlib.sha256()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                hasher.update(chunk)
        return hasher.hexdigest()
    
    def _clean_cache(self):
        """Clean up old versions from cache according to retention policy"""
        cutoff_date = datetime.now() - timedelta(days=self.max_age_days)
        
        for file_name, file_entry in self.cache_index["files"].items():
            versions = file_entry["versions"]
            
            # Keep newest versions based on max_versions
            if len(versions) > self.max_versions:
                for old_version in versions[self.max_versions:]:
                    self._remove_cache_version(old_version)
                file_entry["versions"] = versions[:self.max_versions]
            
            # Remove versions older than max_age_days
            current_versions = file_entry["versions"].copy()
            for version in current_versions:
                try:
                    version_date = datetime.fromisoformat(version["timestamp"])
                    if version_date < cutoff_date:
                        self._remove_cache_version(version)
                        file_entry["versions"].remove(version)
                except (ValueError, KeyError):
                    # Skip versions with invalid timestamps
                    pass
    
    def _remove_cache_version(self, version):
        """Remove a specific version from cache"""
        try:
            # Remove main file
            if os.path.exists(version["path"]):
                os.remove(version["path"])
            
            # Remove optimized versions
            for opt_type, opt_path in version.get("optimized", {}).items():
                if os.path.exists(opt_path):
                    os.remove(opt_path)
            
            logger.info(f"Removed cache version: {version['version_id']}")
        except (IOError, KeyError) as e:
            logger.error(f"Failed to remove cache version {version.get('version_id')}: {e}")
    
    def get_latest_version(self, file_name: str) -> Optional[Dict]:
        """Get info about the latest version of a file"""
        file_entry = self.cache_index["files"].get(file_name, {"versions": []})
        return file_entry["versions"][0] if file_entry["versions"] else None
    
    def get_version_by_id(self, file_name: str, version_id: str) -> Optional[Dict]:
        """Get a specific version of a file by ID"""
        file_entry = self.cache_index["files"].get(file_name, {"versions": []})
        for version in file_entry["versions"]:
            if version["version_id"] == version_id:
                return version
        return None
    
    def list_available_files(self) -> List[str]:
        """List all available files in cache"""
        return list(self.cache_index["files"].keys())
    
    def list_versions(self, file_name: str) -> List[Dict]:
        """List all available versions for a file"""
        file_entry = self.cache_index["files"].get(file_name, {"versions": []})
        return file_entry["versions"]

def main():
    """Command line interface for cache manager"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Manage cached Auspice JSON files")
    parser.add_argument("--update", action="store_true", help="Update cache with latest files")
    parser.add_argument("--force", action="store_true", help="Force update all files")
    parser.add_argument("--cache-dir", default="cache/auspice", help="Cache directory")
    parser.add_argument("--source-dir", default="results/auspice", help="Source directory")
    parser.add_argument("--list", action="store_true", help="List available files")
    parser.add_argument("--list-versions", metavar="FILE", help="List versions for a specific file")
    
    args = parser.parse_args()
    
    cache_manager = AuspiceJsonCache(
        cache_dir=args.cache_dir,
        source_dir=args.source_dir
    )
    
    if args.update:
        cache_manager.update_cache(force=args.force)
        print("Cache updated successfully")
    
    if args.list:
        files = cache_manager.list_available_files()
        print(f"Available files ({len(files)}):")
        for file_name in files:
            latest = cache_manager.get_latest_version(file_name)
            if latest:
                print(f"  {file_name} (latest: {latest['timestamp']})")
            else:
                print(f"  {file_name}")
    
    if args.list_versions:
        versions = cache_manager.list_versions(args.list_versions)
        print(f"Versions for {args.list_versions} ({len(versions)}):")
        for v in versions:
            print(f"  {v['version_id']} - {v['timestamp']} ({v['size']} bytes)")

if __name__ == "__main__":
    main()