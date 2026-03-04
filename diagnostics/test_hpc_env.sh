#!/bin/bash
# Test script for HPC systems that sets HDF5_USE_FILE_LOCKING before Python starts

# CRITICAL: Set this BEFORE python runs
export HDF5_USE_FILE_LOCKING=FALSE

# Now run the Python script
python test_file_validation.py "$@"
