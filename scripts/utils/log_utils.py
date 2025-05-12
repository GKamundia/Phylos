#!/usr/bin/env python3
"""
Enhanced logging utilities for RVF-Nextstrain project
Provides structured logging with consistent formats and contextual information
"""

import os
import sys
import json
import logging
import time
from datetime import datetime
from typing import Dict, Any, Optional, Union

def setup_logger(
    name: str, 
    log_file: Optional[str] = None, 
    level: str = "INFO",
    json_output: bool = False,
    include_timestamp: bool = True,
    log_to_console: bool = True
) -> logging.Logger:
    """
    Set up a logger with consistent formatting and optional JSON output
    
    Args:
        name: Logger name
        log_file: Path to log file (optional)
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        json_output: Whether to format logs as JSON
        include_timestamp: Whether to include timestamp in log format
        log_to_console: Whether to log to console in addition to file
        
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.getLevelName(level))
    
    # Clear any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Create handlers
    handlers = []
    
    # Add file handler if log_file is specified
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        handlers.append(file_handler)
    
    # Add console handler if requested
    if log_to_console:
        console_handler = logging.StreamHandler(sys.stdout)
        handlers.append(console_handler)
    
    # Set formatter based on output type
    if json_output:
        formatter = logging.Formatter('{"timestamp": "%(asctime)s", "name": "%(name)s", "level": "%(levelname)s", "message": %(message)s}')
    else:
        if include_timestamp:
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        else:
            formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    
    # Add formatter to all handlers
    for handler in handlers:
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    return logger

def log_with_context(
    logger: logging.Logger,
    level: str,
    message: str,
    context: Optional[Dict[str, Any]] = None
) -> None:
    """
    Log a message with additional context in a structured format
    
    Args:
        logger: The logger instance
        level: Log level (INFO, WARNING, ERROR, etc.)
        message: Log message
        context: Additional contextual data to include
    """
    log_level = logging.getLevelName(level)
    
    # Create structured log message
    if context:
        # For JSON logs, pass context as part of the message
        if any(hasattr(h.formatter, '_fmt') and '%(message)s' in h.formatter._fmt for h in logger.handlers):
            log_message = json.dumps({
                "msg": message,
                "context": context
            })
        else:
            # For text logs, append context as a formatted string
            context_str = ", ".join(f"{k}={v}" for k, v in context.items())
            log_message = f"{message} [{context_str}]"
    else:
        log_message = message if isinstance(message, str) else str(message)
    
    # Log at the appropriate level
    logger_method = getattr(logger, level.lower())
    logger_method(log_message)

def log_execution_stats(
    logger: logging.Logger, 
    start_time: float, 
    context: Dict[str, Any] = None,
    status: str = "completed"
) -> None:
    """
    Log execution statistics including runtime
    
    Args:
        logger: The logger instance
        start_time: Start time of execution (from time.time())
        context: Additional contextual data
        status: Execution status (completed, failed, etc.)
    """
    end_time = time.time()
    runtime_seconds = end_time - start_time
    
    # Format runtime for human readability
    hours, remainder = divmod(runtime_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    runtime_formatted = f"{int(hours)}h:{int(minutes)}m:{seconds:.2f}s"
    
    # Prepare execution stats
    stats = {
        "runtime_seconds": runtime_seconds,
        "runtime_formatted": runtime_formatted,
        "status": status,
        "timestamp": datetime.now().isoformat()
    }
    
    # Add context if provided
    if context:
        stats.update(context)
    
    log_with_context(
        logger,
        "INFO" if status == "completed" else "ERROR",
        f"Execution {status} in {runtime_formatted}",
        stats
    )

def get_step_logger(step_name: str, pathogen: str, segment: Optional[str] = None) -> logging.Logger:
    """
    Get a logger for a specific pipeline step with standardized naming
    
    Args:
        step_name: Name of the pipeline step
        pathogen: Pathogen name
        segment: Genome segment (if applicable)
        
    Returns:
        Configured logger instance
    """
    log_name = f"{step_name}_{pathogen}"
    if segment:
        log_name = f"{log_name}_{segment}"
    
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{log_name}.log")
    
    return setup_logger(log_name, log_file=log_file)