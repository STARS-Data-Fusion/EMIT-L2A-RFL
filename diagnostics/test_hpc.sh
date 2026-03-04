#!/bin/bash
# HPC-safe wrapper for test_file_validation.py
# This ensures HDF5_USE_FILE_LOCKING is set at the shell level before Python starts

echo "=========================================="
echo "HPC NetCDF Test Script"
echo "=========================================="
echo ""

# Set the critical environment variable
export HDF5_USE_FILE_LOCKING=FALSE
echo "✓ Set HDF5_USE_FILE_LOCKING=FALSE"

# Check if conda environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "EMITL2ARFL" ]]; then
    echo "⚠ WARNING: EMITL2ARFL conda environment not activated"
    echo "  Current environment: ${CONDA_DEFAULT_ENV:-base}"
    echo ""
    echo "  Please activate the environment first:"
    echo "  $ conda activate EMITL2ARFL"
    echo "  $ ./test_hpc.sh"
    echo ""
    exit 1
fi

echo "✓ EMITL2ARFL environment is active"
echo "✓ Python: $(which python)"
echo ""

# Run the test script
echo "Running test_file_validation.py..."
echo "=========================================="
python test_file_validation.py

exit_code=$?

echo ""
echo "=========================================="
if [ $exit_code -eq 0 ]; then
    echo "✓ Test PASSED"
else
    echo "✗ Test FAILED (exit code: $exit_code)"
    echo ""
    echo "Troubleshooting steps:"
    echo "1. Run the diagnostic: python debug_hpc_environment.py <path_to_file>"
    echo "2. Check file permissions: ls -lh ~/data/*.nc"
    echo "3. Verify file integrity: Check file size is not 0"
fi
echo "=========================================="

exit $exit_code
