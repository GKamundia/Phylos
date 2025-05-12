# filepath: c:\Users\Anarchy\Documents\Data_Science\NextStrain\rvf-nextstrain\scripts\caching\api_server.py
#!/usr/bin/env python3
"""
Lightweight API server for serving cached Auspice JSON files
"""

import os
import json
import gzip
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

from fastapi import FastAPI, Response, Query, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from starlette.background import BackgroundTask

from cache_manager import AuspiceJsonCache

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("auspice_api")

app = FastAPI(
    title="RVF-Nextstrain API",
    description="API for serving RVF-Nextstrain Auspice JSON data",
    version="1.0.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Modify for production
    allow_credentials=True,
    allow_methods=["GET", "OPTIONS"],
    allow_headers=["*"],
)

# Cache configuration
CACHE_DIR = os.environ.get("CACHE_DIR", "cache/auspice")
SOURCE_DIR = os.environ.get("SOURCE_DIR", "results/auspice")
ENABLE_COMPRESSION = os.environ.get("ENABLE_COMPRESSION", "1") == "1"

# Initialize cache manager
cache_manager = AuspiceJsonCache(
    cache_dir=CACHE_DIR,
    source_dir=SOURCE_DIR
)

# Request counter for metrics
request_metrics = {
    "total_requests": 0,
    "cache_hits": 0,
    "cache_misses": 0,
    "errors": 0,
    "start_time": datetime.now().isoformat()
}

def should_compress(size_bytes: int, accept_encoding: str) -> bool:
    """Determine if response should be compressed"""
    # Only compress if client accepts gzip and file is large enough
    return (
        ENABLE_COMPRESSION and
        size_bytes > 10000 and  # Only compress files larger than 10KB
        "gzip" in accept_encoding
    )

def compress_json(data: Dict) -> bytes:
    """Compress JSON data using gzip"""
    json_str = json.dumps(data).encode("utf-8")
    return gzip.compress(json_str)

@app.get("/api/v1/dataset/{dataset_name}")
async def get_dataset(
    dataset_name: str,
    version: Optional[str] = None,
    format: str = "full",
    accept_encoding: str = Query("", alias="Accept-Encoding")
):
    """
    Get an Auspice dataset
    
    Args:
        dataset_name: Name of the dataset (without .json extension)
        version: Specific version ID (optional)
        format: Format of the response (full, minimal, tree_only, metadata_only)
        accept_encoding: Accept-Encoding header
    """
    try:
        request_metrics["total_requests"] += 1
        
        # Ensure dataset name ends with .json
        if not dataset_name.endswith(".json"):
            dataset_name = f"{dataset_name}.json"
        
        # Get requested version or latest
        if version:
            dataset_info = cache_manager.get_version_by_id(dataset_name, version)
        else:
            dataset_info = cache_manager.get_latest_version(dataset_name)
        
        if not dataset_info:
            request_metrics["errors"] += 1
            raise HTTPException(status_code=404, detail=f"Dataset {dataset_name} not found")
        
        # Determine file path based on format
        if format == "full":
            file_path = dataset_info["path"]
        elif format in dataset_info.get("optimized", {}):
            file_path = dataset_info["optimized"][format]
        else:
            request_metrics["errors"] += 1
            raise HTTPException(
                status_code=400, 
                detail=f"Invalid format '{format}'. Must be one of: full, minimal, tree_only, metadata_only"
            )
        
        # Check if file exists
        if not os.path.exists(file_path):
            request_metrics["cache_misses"] += 1
            request_metrics["errors"] += 1
            raise HTTPException(status_code=404, detail=f"Dataset file not found")
        
        # Get file size
        file_size = os.path.getsize(file_path)
        
        # Decide whether to compress
        use_compression = should_compress(file_size, accept_encoding)
        
        request_metrics["cache_hits"] += 1
        
        # If compression is enabled and appropriate, load and compress the JSON
        if use_compression:
            try:
                with open(file_path, 'r') as f:
                    data = json.load(f)
                compressed_data = compress_json(data)
                return Response(
                    content=compressed_data,
                    media_type="application/json",
                    headers={"Content-Encoding": "gzip"}
                )
            except Exception as e:
                logger.error(f"Compression error: {e}")
                # Fall back to uncompressed response
        
        # Return the file directly
        return FileResponse(
            path=file_path, 
            media_type="application/json",
            filename=os.path.basename(file_path)
        )
    
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error serving dataset {dataset_name}: {e}")
        request_metrics["errors"] += 1
        raise HTTPException(status_code=500, detail="Internal server error")

@app.get("/api/v1/available")
async def list_available_datasets():
    """List all available datasets"""
    available_files = cache_manager.list_available_files()
    datasets = []
    
    for file_name in available_files:
        latest = cache_manager.get_latest_version(file_name)
        if latest:
            datasets.append({
                "name": file_name,
                "latest_version": latest["version_id"],
                "updated": latest["timestamp"],
                "size": latest["size"],
                "formats": ["full", "minimal", "tree_only", "metadata_only"]
            })
    
    return {"datasets": datasets}

@app.get("/api/v1/versions/{dataset_name}")
async def list_dataset_versions(dataset_name: str):
    """List all versions of a specific dataset"""
    if not dataset_name.endswith(".json"):
        dataset_name = f"{dataset_name}.json"
    
    versions = cache_manager.list_versions(dataset_name)
    if not versions:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_name} not found")
    
    return {"dataset": dataset_name, "versions": versions}

@app.get("/api/v1/health")
async def health_check():
    """API health check endpoint"""
    return {
        "status": "ok",
        "timestamp": datetime.now().isoformat(),
        "metrics": request_metrics
    }

@app.get("/api/v1/update_cache")
async def update_cache(force: bool = False):
    """
    Update the cache with latest Auspice JSON files
    
    This endpoint should be secured in production
    """
    try:
        cache_manager.update_cache(force=force)
        return {"status": "ok", "message": "Cache updated successfully"}
    except Exception as e:
        logger.error(f"Cache update error: {e}")
        raise HTTPException(status_code=500, detail=f"Cache update failed: {str(e)}")

if __name__ == "__main__":
    import uvicorn
    
    # Update cache on startup
    cache_manager.update_cache()
    
    # Start server
    uvicorn.run(app, host="0.0.0.0", port=8000)