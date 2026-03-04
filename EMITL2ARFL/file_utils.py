"""
Utilities for file handling and validation on HPC/network filesystems.

This module provides helper functions for robust file operations in challenging
environments like HPC clusters with network filesystems (NFS, Lustre, etc.).
"""

import logging
import hashlib
import time
from pathlib import Path
from typing import Union, Optional

logger = logging.getLogger(__name__)


def compute_file_checksum(filepath: Union[str, Path], algorithm: str = 'md5', chunk_size: int = 8192) -> str:
    """
    Compute a checksum for a file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the file
    algorithm : str, optional
        Hash algorithm to use ('md5', 'sha256', etc.). Defaults to 'md5'.
    chunk_size : int, optional
        Size of chunks to read at a time (in bytes). Defaults to 8192.
    
    Returns
    -------
    str
        Hexadecimal digest of the file
    """
    hasher = hashlib.new(algorithm)
    filepath = Path(filepath)
    
    with open(filepath, 'rb') as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    
    return hasher.hexdigest()


def wait_for_file_stability(
    filepath: Union[str, Path],
    check_interval: float = 0.5,
    max_checks: int = 10,
    check_size: bool = True
) -> bool:
    """
    Wait for a file to become stable (stop changing).
    
    This is useful on network filesystems where file writes may be buffered
    and not immediately visible, or where downloads may still be in progress.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the file to check
    check_interval : float, optional
        Seconds to wait between stability checks. Defaults to 0.5.
    max_checks : int, optional
        Maximum number of checks before giving up. Defaults to 10.
    check_size : bool, optional
        If True, check file size for stability. Defaults to True.
    
    Returns
    -------
    bool
        True if file became stable, False if max_checks exceeded
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        logger.warning(f"File does not exist: {filepath}")
        return False
    
    previous_size = None
    stable_count = 0
    
    for i in range(max_checks):
        try:
            current_size = filepath.stat().st_size
            
            if check_size and previous_size is not None:
                if current_size == previous_size:
                    stable_count += 1
                    if stable_count >= 2:  # Stable for 2+ checks
                        logger.debug(f"File stable after {i+1} checks: {filepath}")
                        return True
                else:
                    stable_count = 0
            
            previous_size = current_size
            time.sleep(check_interval)
            
        except Exception as e:
            logger.warning(f"Error checking file stability: {e}")
            return False
    
    logger.warning(f"File may not be stable after {max_checks} checks: {filepath}")
    return False


def safe_file_remove(filepath: Union[str, Path], max_attempts: int = 3, delay: float = 0.5) -> bool:
    """
    Safely remove a file with retry logic for network filesystems.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the file to remove
    max_attempts : int, optional
        Maximum number of removal attempts. Defaults to 3.
    delay : float, optional
        Seconds to wait between attempts. Defaults to 0.5.
    
    Returns
    -------
    bool
        True if file was successfully removed or doesn't exist, False otherwise
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        return True
    
    for attempt in range(max_attempts):
        try:
            filepath.unlink()
            time.sleep(delay)  # Wait for filesystem to sync
            
            if not filepath.exists():
                logger.debug(f"Successfully removed file: {filepath}")
                return True
                
        except Exception as e:
            logger.warning(f"Failed to remove file (attempt {attempt + 1}/{max_attempts}): {filepath}. Error: {e}")
            if attempt < max_attempts - 1:
                time.sleep(delay)
    
    logger.error(f"Could not remove file after {max_attempts} attempts: {filepath}")
    return False


def verify_file_readable(filepath: Union[str, Path], read_size: int = 1024) -> bool:
    """
    Verify that a file can be opened and read.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the file to check
    read_size : int, optional
        Number of bytes to attempt reading. Defaults to 1024.
    
    Returns
    -------
    bool
        True if file is readable, False otherwise
    """
    filepath = Path(filepath)
    
    try:
        with open(filepath, 'rb') as f:
            _ = f.read(read_size)
        return True
    except Exception as e:
        logger.warning(f"File not readable: {filepath}. Error: {e}")
        return False
