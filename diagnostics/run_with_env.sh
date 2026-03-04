#!/bin/bash
# Wrapper script to run Python scripts with the correct conda environment
# This ensures netCDF4 and other dependencies are available

# Usage: ./run_with_env.sh script_name.py [arguments...]
# Example: ./run_with_env.sh test_file_validation.py
# Example: ./run_with_env.sh debug_hpc_environment.py ~/data/file.nc

# Find the conda base directory
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
elif [ -f "$HOME/.bash_profile" ]; then
    source "$HOME/.bash_profile"
elif [ -f "$HOME/.zshrc" ]; then
    source "$HOME/.zshrc"
fi

# Try to find conda
if command -v conda &> /dev/null; then
    echo "Activating EMITL2ARFL conda environment..."
    
    # Initialize conda for bash
    eval "$(conda shell.bash hook)"
    
    # Activate the environment
    conda activate EMITL2ARFL
    
    if [ $? -eq 0 ]; then
        echo "Environment activated successfully"
        echo "Python: $(which python)"
        echo ""
        
        # Run the script with all arguments
        python "$@"
    else
        echo "ERROR: Failed to activate EMITL2ARFL environment"
        echo "Make sure it exists: conda env list"
        exit 1
    fi
else
    echo "ERROR: conda not found"
    echo "Please install conda or activate the environment manually:"
    echo "  conda activate EMITL2ARFL"
    echo "  python $@"
    exit 1
fi
