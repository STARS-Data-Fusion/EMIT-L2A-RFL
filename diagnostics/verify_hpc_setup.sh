#!/bin/bash
# Quick verification script for HPC processing workflows
# Run this before starting large processing jobs to ensure environment is correct

echo "=============================================="
echo "EMIT L2A RFL HPC Environment Verification"
echo "=============================================="
echo ""

# Check 1: Environment variable
echo "✓ Checking HDF5_USE_FILE_LOCKING..."
if [ "$HDF5_USE_FILE_LOCKING" = "FALSE" ]; then
    echo "  ✓ HDF5_USE_FILE_LOCKING=FALSE"
else
    echo "  ✗ HDF5_USE_FILE_LOCKING is not set to FALSE"
    echo "  Current value: ${HDF5_USE_FILE_LOCKING:-<not set>}"
    echo ""
    echo "  FIX: Run these commands:"
    echo "    export HDF5_USE_FILE_LOCKING=FALSE"
    echo "    source ~/.bashrc"
    exit 1
fi
echo ""

# Check 2: Conda environment
echo "✓ Checking conda environment..."
if [ "$CONDA_DEFAULT_ENV" = "EMITL2ARFL" ]; then
    echo "  ✓ EMITL2ARFL environment active"
    echo "  Python: $(which python)"
else
    echo "  ✗ EMITL2ARFL environment not active"
    echo "  Current: ${CONDA_DEFAULT_ENV:-none}"
    echo ""
    echo "  FIX: conda activate EMITL2ARFL"
    exit 1
fi
echo ""

# Check 3: Python imports
echo "✓ Checking Python imports..."
python -c "
import sys
try:
    import netCDF4
    print(f'  ✓ netCDF4 {netCDF4.__version__}')
except ImportError as e:
    print(f'  ✗ Cannot import netCDF4: {e}')
    sys.exit(1)

try:
    import EMITL2ARFL
    print(f'  ✓ EMITL2ARFL {EMITL2ARFL.__version__}')
except ImportError as e:
    print(f'  ✗ Cannot import EMITL2ARFL: {e}')
    sys.exit(1)

try:
    import xarray
    print(f'  ✓ xarray {xarray.__version__}')
except ImportError:
    print('  ⚠ xarray not available (optional)')

try:
    import earthaccess
    print(f'  ✓ earthaccess available')
except ImportError:
    print('  ⚠ earthaccess not available (optional)')
" || exit 1

echo ""

# Check 4: Test file access (if file exists)
if [ -n "$1" ] && [ -f "$1" ]; then
    echo "✓ Testing NetCDF file access: $1"
    python -c "
import sys
from EMITL2ARFL import validate_NetCDF_file
try:
    validate_NetCDF_file('$1')
    print('  ✓ File validation successful')
except Exception as e:
    print(f'  ✗ File validation failed: {e}')
    sys.exit(1)
" || exit 1
    echo ""
elif [ -n "$1" ]; then
    echo "⚠ Warning: Test file not found: $1"
    echo ""
fi

# Summary
echo "=============================================="
echo "✓ ALL CHECKS PASSED"
echo "=============================================="
echo ""
echo "Your HPC environment is properly configured."
echo "You can now run your processing scripts:"
echo "  python generate_kings_canyon_timeseries.py"
echo "  python generate_EMIT_L2A_RFL_timeseries.py"
echo ""
echo "For batch jobs, ensure your job script includes:"
echo "  export HDF5_USE_FILE_LOCKING=FALSE"
echo "=============================================="
