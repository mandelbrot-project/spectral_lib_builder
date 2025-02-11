#!/bin/bash
#SBATCH --partition=public-cpu
#SBATCH --time=00:30
#SBATCH --mail-user=adrutz@ethz.ch
#SBATCH --mail-type=ALL

# Create unique working directory in $TMPDIR for this job
WORK_DIR="${TMPDIR}/cfm_job_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${WORK_DIR}"
mkdir -p "${WORK_DIR}/testout"

# Get input file
FILE_INDEX=$(printf "%05d" ${SLURM_ARRAY_TASK_ID:-1})
FILE="test/test-${FILE_INDEX}.txt"

# Copy necessary files to scratch space
cp "${FILE}" "${WORK_DIR}/"
cp -r cfm-4 "${WORK_DIR}/"

# Process the file in scratch space
cd "${WORK_DIR}"

# Run CFM-ID with output going to scratch space
singularity run cfm-4/cfmid_latest.sif cfm-predict \
    "test-${FILE_INDEX}.txt" \
    0.001 \
    /trained_models_cfmid4.0/[M+H]+/param_output.log \
    /trained_models_cfmid4.0/[M+H]+/param_config.txt \
    1 \
    testout \
    0 \
    0

# Create final output directory if it doesn't exist
mkdir -p "$SLURM_SUBMIT_DIR/results/testout"

# Copy results from scratch to permanent storage
cp -r testout/* "$SLURM_SUBMIT_DIR/results/testout/"