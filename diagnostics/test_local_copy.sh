#!/bin/bash
# Test if the file works when copied to local scratch storage

export HDF5_USE_FILE_LOCKING=FALSE

if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_netcdf_file>"
    exit 1
fi

SOURCE_FILE="$1"
TEMP_FILE="/tmp/test_emit_$(basename "$SOURCE_FILE")"

echo "Testing file on network storage..."
echo "Source: $SOURCE_FILE"
echo ""

# Test 1: Try to read from network storage
echo "=== Test 1: Read from network storage ==="
python test_file_validation.py "$SOURCE_FILE"
NETWORK_RESULT=$?
echo ""

# Test 2: Copy to local /tmp and try again
echo "=== Test 2: Copy to local scratch and test ==="
echo "Copying to: $TEMP_FILE"
cp "$SOURCE_FILE" "$TEMP_FILE"
echo "Testing local copy..."
python test_file_validation.py "$TEMP_FILE"
LOCAL_RESULT=$?
echo ""

# Cleanup
rm -f "$TEMP_FILE"

# Summary
echo "=== SUMMARY ==="
echo "Network storage result: $([ $NETWORK_RESULT -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "Local scratch result:   $([ $LOCAL_RESULT -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo ""

if [ $NETWORK_RESULT -ne 0 ] && [ $LOCAL_RESULT -eq 0 ]; then
    echo "⚠️  DIAGNOSIS: Network filesystem caching/metadata issue"
    echo ""
    echo "The file is NOT corrupted - it's a network filesystem problem."
    echo ""
    echo "SOLUTIONS:"
    echo "1. Download directly to local scratch first:"
    echo "   python download_to_scratch.py"
    echo ""
    echo "2. Add to your ~/.config/fish/config.fish:"
    echo "   set -Ux HDF5_USE_FILE_LOCKING FALSE"
    echo ""
    echo "3. For existing files, copy to local scratch before processing:"
    echo "   cp <network_file> /tmp/ && process_from_tmp"
    echo ""
    echo "4. Check filesystem type with: df -T $SOURCE_FILE"
elif [ $NETWORK_RESULT -eq 0 ]; then
    echo "✓ File works on network storage - no issues detected"
else
    echo "✗ File appears genuinely corrupted - try re-downloading"
fi
