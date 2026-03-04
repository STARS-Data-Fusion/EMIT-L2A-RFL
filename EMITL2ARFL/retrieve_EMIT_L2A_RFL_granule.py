import posixpath
import logging
import time
from os import remove, sync
from os.path import join, expanduser, abspath, exists
from typing import List, Optional

import earthaccess

from .constants import *
from .EMITL2ARFLGranule import EMITL2ARFLGranule
from .find_EMIT_L2A_RFL_granule import find_EMIT_L2A_RFL_granule
from .validate_NetCDF_file import validate_NetCDF_file
from .exceptions import NetCDFValidationError
from .file_utils import safe_file_remove, wait_for_file_stability

logger = logging.getLogger(__name__)

def retrieve_EMIT_L2A_RFL_granule(
        remote_granule: earthaccess.search.DataGranule = None,
        orbit: int = None,
        scene: int = None, 
        download_directory: str = DOWNLOAD_DIRECTORY,
        max_retries: int = 3,
        retry_delay: float = 2.0,
        skip_validation: bool = False,
        threads: int = 1,
        use_wget: bool = False) -> EMITL2ARFLGranule:
    """
    Retrieve an EMIT L2A Reflectance granule with resilient error handling and retry logic.

    This function retrieves an EMIT L2A Reflectance granule based on the provided granule, orbit, and scene.
    If the granule is not provided, it searches for the granule using the orbit and scene parameters.
    The granule is then downloaded to the specified directory and validated. If files are missing or corrupted,
    the function will attempt to re-download them up to max_retries times.

    Args:
        remote_granule (earthaccess.search.DataGranule, optional): The granule to retrieve. Defaults to None.
        orbit (int, optional): The orbit number to search for the granule. Defaults to None.
        scene (int, optional): The scene number to search for the granule. Defaults to None.
        download_directory (str, optional): The directory to download the granule files to. Defaults to DOWNLOAD_DIRECTORY.
        max_retries (int, optional): Maximum number of retry attempts for downloading corrupted files. Defaults to 3.
        retry_delay (float, optional): Seconds to wait between retry attempts. Useful for HPC environments. Defaults to 2.0.
        skip_validation (bool, optional): If True, skip NetCDF validation (use with caution). Defaults to False.
            Only use this if you're experiencing persistent corruption and want to attempt processing anyway.
        threads (int, optional): Number of parallel download threads. Defaults to 1 (single-threaded) for maximum
            reliability on HPC/network filesystems. Set to 8 for faster downloads on local systems with stable storage.
        use_wget (bool, optional): If True, use wget command-line tool for downloading. If False (default), use 
            earthaccess.download(). Note: wget approach is currently not implemented for NASA Earthdata authentication.
            This parameter is reserved for future use.

    Returns:
        EMITL2ARFLGranule: The retrieved EMIT L2A Reflectance granule wrapped in an EMITL2ARFLGranule object.

    Raises:
        ValueError: If no granule is found for the provided orbit and scene, or if the provided granule is not an EMIT L2A Reflectance collection 1 granule.
        FileNotFoundError: If required files cannot be downloaded after max_retries attempts.
        RuntimeError: If use_wget=True but wget is not available on the system.
    """
    if remote_granule is None and orbit is not None and scene is not None:
        remote_granule = find_EMIT_L2A_RFL_granule(granule=remote_granule, orbit=orbit, scene=scene)
    
    if remote_granule is None:
        raise ValueError("either granule or orbit and scene must be provided")

    # Parse granule ID from the first data link
    granule_ID = posixpath.splitext(posixpath.basename(remote_granule.data_links()[0]))[0]

    # Validate that this is an EMIT L2A Reflectance collection 1 granule
    if not granule_ID.startswith("EMIT_L2A_RFL_001_"):
        raise ValueError("The provided granule is not an EMIT L2A Reflectance collection 1 granule.")

    # Set up the granule directory
    directory = join(download_directory, granule_ID)
    abs_directory = abspath(expanduser(directory))
    
    # Get expected filenames from the remote granule
    base_filenames = [posixpath.basename(URL) for URL in remote_granule.data_links()]
    local_files = [join(abs_directory, fname) for fname in base_filenames]
    
    # Log the URLs that will be checked/downloaded
    logger.info(f"Granule URLs:")
    for i, url in enumerate(remote_granule.data_links(), 1):
        logger.info(f"  [{i}/{len(remote_granule.data_links())}] {url}")
    
    # Helper function to identify file types
    def _identify_files(files: List[str]) -> tuple[Optional[str], Optional[str], Optional[str]]:
        """Identify RFL, MASK, and RFLUNCERT files from a list of filenames."""
        reflectance = next((f for f in files if '_RFL_' in f and '_RFLUNCERT_' not in f), None)
        mask = next((f for f in files if '_MASK_' in f), None)
        uncertainty = next((f for f in files if '_RFLUNCERT_' in f), None)
        return reflectance, mask, uncertainty
    
    # Helper function to download specific files
    def _download_files(urls: List[str], retry_attempt: int = 0) -> bool:
        """Download files using earthaccess (with optional wget support for future use)."""
        try:
            if use_wget:
                # TODO: Implement wget support with NASA Earthdata authentication
                # For now, fall back to earthaccess
                logger.warning("wget download not yet implemented for NASA Earthdata - using earthaccess")
            
            # Use earthaccess for downloading (handles NASA authentication)
            actual_threads = threads
            logger.info(f"Downloading with earthaccess (attempt {retry_attempt + 1}/{max_retries}, threads={actual_threads})...")
            logger.info(f"Download directory: {abs_directory}")
            for i, url in enumerate(urls, 1):
                logger.info(f"  [{i}/{len(urls)}] {url}")
            
            earthaccess.download(urls, local_path=abs_directory, threads=actual_threads)
            return True
        except Exception as e:
            logger.error(f"Download failed: {e}")
            return False
    
    # Identify expected file paths
    reflectance_filename, mask_filename, uncertainty_filename = _identify_files(local_files)
    
    if not all([reflectance_filename, mask_filename, uncertainty_filename]):
        raise ValueError('Could not identify all required file types (RFL, MASK, RFLUNCERT) from granule data links.')
    
    # Track which files need to be (re)downloaded
    files_to_download = []
    file_info = {
        'reflectance': (reflectance_filename, 'Reflectance'),
        'mask': (mask_filename, 'Mask'),
        'uncertainty': (uncertainty_filename, 'Uncertainty')
    }
    
    # Initial validation check (unless validation is skipped)
    if not skip_validation:
        for filepath, file_type in file_info.values():
            try:
                validate_NetCDF_file(filepath, file_type=file_type)
                logger.info(f"Cached file validated successfully: {filepath}")
            except NetCDFValidationError as e:
                logger.warning(f"Cached file validation failed: {e}")
                # Remove corrupted cached file immediately to force re-download
                if exists(filepath):
                    logger.info(f"Removing corrupted cached file: {filepath}")
                    if safe_file_remove(filepath, max_attempts=3):
                        logger.info(f"Corrupted file removed successfully")
                    else:
                        logger.error(f"Could not remove corrupted file - will attempt download anyway")
                files_to_download.append(filepath)
    else:
        logger.warning("Validation is SKIPPED - files may be corrupted!")
        # Check if files exist but don't validate them
        for filepath, file_type in file_info.values():
            if not exists(filepath):
                files_to_download.append(filepath)
    
    # Retry loop for downloading/repairing files
    retry_count = 0
    while files_to_download and retry_count < max_retries:
        # Wait before retry (except on first attempt)
        if retry_count > 0:
            wait_time = retry_delay * (2 ** (retry_count - 1))  # Exponential backoff
            logger.info(f"Waiting {wait_time:.1f} seconds before retry {retry_count + 1}...")
            time.sleep(wait_time)
        
        # Force filesystem sync (important for HPC/network filesystems)
        try:
            sync()
        except Exception:
            pass  # sync() may not be available on all systems
        
        # Download all files (earthaccess will skip existing valid files)
        if not _download_files(remote_granule.data_links(), retry_count):
            retry_count += 1
            continue
        
        # Wait for files to stabilize (important for network filesystems)
        logger.info("Waiting for downloaded files to stabilize...")
        for filepath in [reflectance_filename, mask_filename, uncertainty_filename]:
            if exists(filepath):
                wait_for_file_stability(filepath, check_interval=0.5, max_checks=6)
        
        # Force filesystem sync after download
        try:
            sync()
            time.sleep(0.5)  # Brief pause to ensure file writes complete
        except Exception:
            pass
        
        # Re-validate files (unless validation is skipped)
        files_to_download = []
        if not skip_validation:
            for filepath, file_type in file_info.values():
                try:
                    validate_NetCDF_file(filepath, file_type=file_type)
                    logger.info(f"Downloaded file validated successfully: {filepath}")
                except NetCDFValidationError as e:
                    logger.warning(f"Validation failed after download attempt: {e}")
                    # File still corrupted after download - remove it for next retry
                    if exists(filepath):
                        logger.info(f"Removing still-corrupted file for retry: {filepath}")
                        safe_file_remove(filepath, max_attempts=3)
                    files_to_download.append(filepath)
        
        if not files_to_download:
            logger.info("All files successfully downloaded and validated.")
            break
        
        retry_count += 1
    
    # Final validation - raise error if any files are still invalid (unless validation is skipped)
    if files_to_download and not skip_validation:
        error_msg = (
            f"Failed to download valid files after {max_retries} attempts. "
            f"Missing or corrupted files: {files_to_download}\n\n"
            f"TROUBLESHOOTING STEPS:\n"
            f"1. Files are consistently corrupted - this suggests an HPC/network filesystem issue\n"
            f"2. Try downloading to local scratch space instead of network storage:\n"
            f"   download_directory='/tmp/emit_download'\n"
            f"3. Increase retries and delays:\n"
            f"   max_retries=10, retry_delay=10.0\n"
            f"4. Try single-threaded downloads by setting threads=1 in earthaccess.download()\n"
            f"5. As a last resort, you can skip validation (files may still be unusable):\n"
            f"   retrieve_EMIT_L2A_RFL_granule(..., skip_validation=True)\n"
            f"6. Contact your HPC support team - this may indicate storage system issues"
        )
        raise FileNotFoundError(error_msg)
    
    # Create and return the granule object
    local_granule = EMITL2ARFLGranule(
        reflectance_filename=reflectance_filename,
        mask_filename=mask_filename,
        uncertainty_filename=uncertainty_filename
    )

    return local_granule