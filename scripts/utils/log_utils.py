#!/usr/bin/env python3
"""
Structured logging utilities for the RVF-Nextstrain pipeline
Provides consistent formatting and logging levels across all scripts
"""

import os
import sys
import json
import logging
from datetime import datetime
from pathlib import Path

# Default log format with color support for terminal and structured format for files
DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
JSON_FORMAT = {
    "timestamp": "%(asctime)s",
    "level": "%(levelname)s",
    "name": "%(name)s",
    "module": "%(module)s",
    "function": "%(funcName)s",
    "line": "%(lineno)d",
    "message": "%(message)s"
}

# Define log levels for different types of events
LOG_LEVELS = {
    "DEBUG": logging.DEBUG,      # Detailed diagnostic information
    "INFO": logging.INFO,        # Confirmation that things are working
    "WARNING": logging.WARNING,  # Indication that something unexpected happened
    "ERROR": logging.ERROR,      # The software has not been able to perform a task
    "CRITICAL": logging.CRITICAL # A serious error that may prevent the program from continuing
}

class StructuredFormatter(logging.Formatter):
    """Custom formatter for structured logs that can output both text and JSON"""
    
    def __init__(self, fmt=None, json_fmt=None, datefmt=None, json_output=False):
        super().__init__(fmt, datefmt)
        self.json_fmt = json_fmt or JSON_FORMAT
        self.json_output = json_output
    
    def format(self, record):
        if self.json_output:
            log_data = {key: value % vars(record) 
                        for key, value in self.json_fmt.items()
                        if hasattr(record, key.split("(")[0].strip("%"))}
            
            # Add any extra attributes
            if hasattr(record, "extra"):
                log_data.update(record.extra)
            
            # Add exception info if present
            if record.exc_info:
                log_data["exception"] = {
                    "type": record.exc_info[0].__name__,
                    "message": str(record.exc_info[1]),
                    "traceback": self.formatException(record.exc_info)
                }
                
            return json.dumps(log_data)
        else:
            return super().format(record)

def setup_logger(name, 
                log_file=None, 
                level="INFO", 
                json_output=False, 
                console_output=True,
                log_dir="logs"):
    """
    Sets up a logger with both file and console handlers
    
    Args:
        name (str): Name of the logger
        log_file (str, optional): Path to log file. If None, no file handler is created.
        level (str, optional): Log level. Defaults to "INFO".
        json_output (bool, optional): Whether to output logs in JSON format. Defaults to False.
        console_output (bool, optional): Whether to output logs to console. Defaults to True.
        log_dir (str, optional): Directory to store log files. Defaults to "logs".
        
    Returns:
        logging.Logger: Configured logger instance
    """
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(LOG_LEVELS.get(level.upper(), logging.INFO))
    
    # Remove existing handlers to avoid duplicates
    logger.handlers = []
    
    # Create formatters
    text_formatter = StructuredFormatter(DEFAULT_FORMAT, json_output=False)
    json_formatter = StructuredFormatter(json_fmt=JSON_FORMAT, json_output=True)
    
    # Add console handler if requested
    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(text_formatter)
        logger.addHandler(console_handler)
    
    # Add file handler if log_file is provided
    if log_file:
        # Ensure log directory exists
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        
        # Use specified log file or create one based on script name
        log_path = log_file
        if not log_path.startswith(log_dir):
            log_path = os.path.join(log_dir, log_file)
            
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(json_formatter if json_output else text_formatter)
        logger.addHandler(file_handler)
    
    return logger

def log_execution_stats(logger, start_time, metadata=None, status="completed"):
    """
    Log execution statistics including runtime
    
    Args:
        logger (logging.Logger): Logger instance
        start_time (float): Start time from time.time()
        metadata (dict, optional): Additional metadata to log
        status (str, optional): Execution status. Defaults to "completed".
    """
    import time
    
    end_time = time.time()
    runtime_seconds = end_time - start_time
    hours, rem = divmod(runtime_seconds, 3600)
    minutes, seconds = divmod(rem, 60)
    
    stats = {
        "status": status,
        "runtime": {
            "seconds": runtime_seconds,
            "formatted": f"{int(hours):02d}:{int(minutes):02d}:{seconds:.2f}"
        }
    }
    
    if metadata:
        stats.update(metadata)
    
    if status == "completed":
        logger.info(f"Execution completed in {stats['runtime']['formatted']}", extra={"stats": stats})
    else:
        logger.warning(f"Execution {status} after {stats['runtime']['formatted']}", extra={"stats": stats})

def log_with_context(logger, level, message, context=None):
    """
    Log a message with additional context
    
    Args:
        logger (logging.Logger): Logger instance
        level (str): Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        message (str): Log message
        context (dict, optional): Additional context to include
    """
    log_level = LOG_LEVELS.get(level.upper(), logging.INFO)
    extra = {"extra": context} if context else {}
    
    logger.log(log_level, message, extra=extra)