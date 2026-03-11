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

# Navigate to project directory
cd $HOME/data/EMIT-L2A-RFL

# Run the dataset generation script
echo "Starting EMIT dataset generation for 2024 February though July..."
python -u generate_kings_canyon_timeseries_2024_feb_july.py

echo "EMIT dataset generation for 2024 completed with exit code: $?"
