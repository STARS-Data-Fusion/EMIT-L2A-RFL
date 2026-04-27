#!/bin/bash

#SBATCH --job-name=hls-kings-2024
#SBATCH --output=logs/hls_kings_2024_%j.log
#SBATCH --error=logs/hls_kings_2024_%j.err
#SBATCH --open-mode=append
#SBATCH --time=48:00:00
#SBATCH --mem=2GB
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${USER}@example.com

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate EMITL2ARFL

# Critical for HPC/network filesystems when reading NetCDF/HDF5 files
export HDF5_USE_FILE_LOCKING=FALSE

# Navigate to project directory
cd $HOME/data/EMIT-L2A-RFL

# Run the dataset generation script
echo "Starting EMIT dataset generation for 2022 August through 2025 November..."
python -u generate_kings_canyon_timeseries.py

echo "EMIT dataset generation for 2024 completed with exit code: $?"
