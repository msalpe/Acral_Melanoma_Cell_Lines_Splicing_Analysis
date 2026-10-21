#!/bin/bash

# ===============================================
# SLURM Job Array Script for Vast-Tools Align
# ===============================================

##################
# SLURM settings #
##################

#SBATCH --output=logs/VastAlign_%A_%a.out
#SBATCH --error=logs/VastAlign_%A_%a.err
#SBATCH --qos=long
#SBATCH --job-name=VastAlign
#SBATCH --time=01-00:00:00    # max por long qos
#SBATCH --cpus-per-task=8     # CPUs per task
#SBATCH --mem=32G             # Memory allocation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-11 #12 samples

set -euo pipefail

# Define directories
PROJECT_DIR="${PROJECT_DIR:?ERROR: PROJECT_DIR is not defined}"
FASTQ_DIR="${PROJECT_DIR}/FASTQ"
OUTPUT_DIR="${PROJECT_DIR}/vast_out"

# suffixes for FASTQ files
R1_SUFFIX="_1.fastq.gz"
R2_SUFFIX="_2.fastq.gz"

##########################################################################################

# Load necessary modules
module load R/4.4.0-gfbf-2023b

# Ensure output directory exists
if [[ ! -d "$OUTPUT_DIR" ]]; then
  echo "Creating output directory: $OUTPUT_DIR"
  mkdir -p "$OUTPUT_DIR"
fi

# Change working directory to output
cd "$OUTPUT_DIR" || { echo "ERROR: Cannot change wd to $OUTPUT_DIR"; exit 1; }

# Ensure FASTQ directory exists
if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "ERROR: FASTQ directory does not exist: $FASTQ_DIR"
  exit 1
fi

# Initialize an array to store sample names and paths
SAMPLES=()

# Recursively search for all FASTQ files ending in "R1 suffix" within FASTQ_DIR,
# including its subdirectories. The -print0 option ensures that filenames with spaces 
# are handled correctly. These will later be paired with their corresponding "R2 suffix" files.
while IFS= read -r -d '' FASTQ1; do
  SAMPLE_NAME=$(basename "$FASTQ1" "$R1_SUFFIX")
  SAMPLE_DIR=$(dirname "$FASTQ1")
  FASTQ2="${SAMPLE_DIR}/${SAMPLE_NAME}${R2_SUFFIX}"

  # Ensure both files exist before adding to the sample list
  if [[ -f "$FASTQ2" ]]; then
    SAMPLES+=("$SAMPLE_NAME|$SAMPLE_DIR")
  else
    echo "WARNING: No matching paired file for sample $SAMPLE_NAME in $SAMPLE_DIR"
  fi
done < <(find "$FASTQ_DIR" -type f -name "*${R1_SUFFIX}" -print0)
# Get the total number of available samples
TOTAL_SAMPLES=${#SAMPLES[@]}

# Ensure there are samples to process
if [[ $TOTAL_SAMPLES -eq 0 ]]; then
  echo "ERROR: No FASTQ samples found in $FASTQ_DIR"
  exit 1
fi

# Ensure SLURM_ARRAY_TASK_ID is within the sample range
if [[ $SLURM_ARRAY_TASK_ID -ge $TOTAL_SAMPLES ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds available samples ($TOTAL_SAMPLES)."
  exit 1
fi

# Select sample for the current job
IFS='|' read -r SAMPLE_NAME SAMPLE_DIR <<< "${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# Define absolute paths for FASTQ files
FASTQ1="${SAMPLE_DIR}/${SAMPLE_NAME}${R1_SUFFIX}"

FASTQ2="${SAMPLE_DIR}/${SAMPLE_NAME}${R2_SUFFIX}"


# Verify that both FASTQ files exist before running Vast-tools align
if [[ -f "$FASTQ1" && -f "$FASTQ2" ]]; then
  echo "Processing sample: $SAMPLE_NAME in $SAMPLE_DIR"
  vast-tools align "$FASTQ1" "$FASTQ2" -sp Hs2 -o "$OUTPUT_DIR" --IR_version 2 --expr -c 8 -n "$SAMPLE_NAME"
else
  echo "ERROR: Missing FASTQ files for sample: $SAMPLE_NAME"
  exit 1
fi
