#!/bin/bash

# ================================================================================
# SLURM Job ARRAY Script for fastqc
# Sequencing FASTQ QC
# ================================================================================

#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err
#SBATCH --qos=vshort
#SBATCH --job-name=FASTQC
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-23 #24 files (12 samples)

set -euo pipefail

# === Define paths ===
PROJECT_DIR="${PROJECT_DIR:?ERROR: PROJECT_DIR is not defined}"
OUT_DIR="${PROJECT_DIR}/fastqc"
FASTQ_LIST="${PROJECT_DIR}/code/fastq_list.txt"

####################################################################################
# Load FastQC
module load FastQC/0.12.1-Java-11

mkdir -p "$OUT_DIR"

# fastqc for each file in the array
FASTQ_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$FASTQ_LIST")

echo "[$SLURM_ARRAY_TASK_ID] Processing: $FASTQ_FILE"

fastqc -o "$OUT_DIR" -t 2 "$FASTQ_FILE"

echo "[$(date)] [$SLURM_ARRAY_TASK_ID] Finished processing."
