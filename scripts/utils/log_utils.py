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
from logging.handlers import RotatingFileHandler

class JsonFormatter(logging.Formatter):
    """
    Custom JSON formatter for structured logging.
    """
    def format(self, record: logging.LogRecord) -> str:
        log_record = {
            "time": self.formatTime(record, self.datefmt),
            "level": record.levelname,
            "name": record.name,
            "message": record.getMessage(),
        }
        if record.exc_info:
            log_record["exception"] = self.formatException(record.exc_info)
        return json.dumps(log_record)

def setup_logger(
    name: str,
    log_file: str = None,
    level: str = "INFO",
    console_output: bool = True,
    json_output: bool = False,
    max_bytes: int = 10*1024*1024, # 10 MB
    backup_count: int = 5
) -> logging.Logger:
    """
    Set up a logger with configurable handlers and formatting.
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.propagate = False # Prevent duplicate logs if root logger is configured

    # Clear existing handlers for this logger to avoid duplication if called multiple times
    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = JsonFormatter() if json_output else logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    if log_file:
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir): # Check if log_dir is not empty before creating
            os.makedirs(log_dir, exist_ok=True)
        
        # Use RotatingFileHandler for better log management
        file_handler = RotatingFileHandler(
            log_file, 
            maxBytes=max_bytes, 
            backupCount=backup_count
        )
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
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