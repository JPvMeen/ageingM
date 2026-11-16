#!/bin/bash
#SBATCH --job-name=influenza_array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=10-00:00:00
#SBATCH --chdir=/hpc/dla_lti/jvanmeenen/ageingM
#SBATCH --output=logs/influenza_array_%A_%a.log

set -euo pipefail
: "${INFLUENZA_COVERAGE:?Set INFLUENZA_COVERAGE to 30 or 80}"
: "${SLURM_ARRAY_JOB_ID:?Run as a Slurm job array}"
: "${SLURM_ARRAY_TASK_ID:?Missing Slurm array task index}"

export OPENBLAS_NUM_THREADS=1
export JULIA_CPU_TARGET=generic
export JULIA_NUM_PRECOMPILE_TASKS=1

echo "Start: $(date); host: $(hostname); coverage: ${INFLUENZA_COVERAGE}%; forcing task: ${SLURM_ARRAY_TASK_ID}; batch: ${SLURM_ARRAY_JOB_ID}"
srun "$HOME/.juliaup/bin/julia" \
    --startup-file=no \
    --project=. \
    --threads="${SLURM_CPUS_PER_TASK}" \
    scripts/runs/ageing_influenza_array.jl
echo "Array task finished at $(date)"
