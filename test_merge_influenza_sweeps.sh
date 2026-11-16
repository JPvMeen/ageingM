#!/bin/bash
#SBATCH --job-name=influenza_merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/hpc/dla_lti/jvanmeenen/ageingM
#SBATCH --output=logs/influenza_merge_%j.log

set -euo pipefail
: "${INFLUENZA_COVERAGE:?Set INFLUENZA_COVERAGE to 30 or 80}"
: "${INFLUENZA_BATCH_ID:?Set INFLUENZA_BATCH_ID to Slurm array job ID}"
export OPENBLAS_NUM_THREADS=1
export JULIA_CPU_TARGET=generic
export JULIA_NUM_PRECOMPILE_TASKS=1
echo "Merging coverage ${INFLUENZA_COVERAGE}% batch ${INFLUENZA_BATCH_ID} at $(date)"
srun "$HOME/.juliaup/bin/julia" \
    --startup-file=no \
    --project=. \
    --threads=1 \
    scripts/runs/merge_influenza_sweeps.jl "${INFLUENZA_COVERAGE}" "${INFLUENZA_BATCH_ID}"
echo "Merge complete at $(date)"
