#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=logs/VastCombine_%j.out
#SBATCH --error=logs/VastCombine.%j.err

# queue
#SBATCH --qos=shorter
#SBATCH --job-name=VastCombine
#SBATCH --time=3:00:00             # max time for short queue
#SBATCH --mem=24G                  # mem 
#SBATCH --cpus-per-task=8          # cpus
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail

# Define directories
PROJECT_DIR="${PROJECT_DIR:?ERROR: PROJECT_DIR is not defined}"
OUTPUT_DIR="${PROJECT_DIR}/vast_out"

# Change working directory
cd "$OUTPUT_DIR" || { echo "ERROR: Cannot access $OUTPUT_DIR"; exit 1; }

# load modules
module load R/4.4.0-gfbf-2023b

# Running Vast-tools combine
vast-tools merge --groups merge_config.txt -sp Hs2 --expr -o "$OUTPUT_DIR"   
vast-tools combine -o "$OUTPUT_DIR" -sp Hs2 -C --add_version
