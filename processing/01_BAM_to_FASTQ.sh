#!/bin/bash

# ================================================================================
# SLURM Job array for converting BAM files to FASTQ using samtools
# ================================================================================

#SBATCH --output=logs/samtools_%A_%a.out
#SBATCH --error=logs/samtools_%A_%a.err
#SBATCH --qos=shorter
#SBATCH --job-name=samtools
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-11  #12samples

set -euo pipefail

# === Define paths ===
PROJECT_DIR="${PROJECT_DIR:?ERROR: PROJECT_DIR is not defined}"
BAM_DIR="${PROJECT_DIR}/BAMs/DATA"
OUT_DIR="${PROJECT_DIR}/FASTQ"

# === Create output directories if they don't exist ===
mkdir -p "$OUT_DIR"


# === Load SAMTools modules ===
module load SAMtools/1.20-GCC-13.2.0 

########################################################################################

# === Map Array ID to File ===
# Generate a list of BAM files and select the one corresponding to the current task ID
mapfile -t BAM_FILES < <(find "$BAM_DIR" -type f -name "*.bam" | sort)
INPUT_BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# Extract base name (remove path and .bam extension)
SAMPLE_NAME=$(basename "$INPUT_BAM" .bam)

# Define output filenames
R1_OUT="${OUT_DIR}/${SAMPLE_NAME}_1.fastq.gz"
R2_OUT="${OUT_DIR}/${SAMPLE_NAME}_2.fastq.gz"

echo "----------------------------------------------------------------"
echo "Processing sample: $SAMPLE_NAME"
echo "Input BAM: $INPUT_BAM"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Date: $(date)"
echo "----------------------------------------------------------------"

# === Run Samtools conversion ===
# 1. samtools collate: Groups reads by name (required for coordinate-sorted BAMs) 
#    to ensure R1 and R2 stay paired.
# 2. samtools fastq: Converts the streamed BAM data into compressed FASTQ files.
# Using pipes and the number of threads defined in the SLURM headers.

samtools collate -u -O -@ $SLURM_CPUS_PER_TASK "$INPUT_BAM" | \
samtools fastq -@ $SLURM_CPUS_PER_TASK \
    -1 "$R1_OUT" \
    -2 "$R2_OUT" \
    -0 /dev/null -s /dev/null -n

echo "----------------------------------------------------------------"
echo "Conversion finished for $SAMPLE_NAME at $(date)"
echo "----------------------------------------------------------------"
