#!/usr/bin/env python3
"""
Retry utilities for handling transient failures
"""

import time
import logging
import functools
import random
from typing import Callable, Any, Optional, Type, Union, List

logger = logging.getLogger(__name__)

def retry(
    max_tries: int = 3,
    exceptions: Union[Type[Exception], List[Type[Exception]]] = Exception,
    delay_seconds: float = 1.0,
    backoff: float = 2.0,
    jitter: float = 0.1,
    logger: Optional[logging.Logger] = None
) -> Callable:
    """
    Retry decorator with exponential backoff for handling transient errors
    
    Args:
        max_tries: Maximum number of attempts
        exceptions: Exception(s) to catch and retry on
        delay_seconds: Initial delay between retries in seconds
        backoff: Multiplier for delay between retries
        jitter: Random factor to add to delay (0.1 = ±10%)
        logger: Logger to use (defaults to module logger)
    
    Returns:
        Decorator function
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            log = logger or logging.getLogger(func.__module__)
            tries = 0
            delay = delay_seconds
            
            while tries < max_tries:
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    tries += 1
                    if tries == max_tries:
                        log.error(f"Failed after {max_tries} attempts: {str(e)}")
                        raise
                    
                    # Calculate delay with jitter
                    jitter_range = delay * jitter
                    actual_delay = delay + random.uniform(-jitter_range, jitter_range)
                    actual_delay = max(0.1, actual_delay)  # Ensure positive delay
                    
                    log.warning(
                        f"Attempt {tries}/{max_tries} failed: {str(e)}. "
                        f"Retrying in {actual_delay:.2f} seconds..."
                    )
                    
                    time.sleep(actual_delay)
                    delay *= backoff
            
            return None  # Should never reach here
        
        return wrapper
    
    return decorator