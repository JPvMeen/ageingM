#!/bin/bash
#SBATCH --job-name=influenza_ts_fva91
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --chdir=/hpc/dla_lti/jvanmeenen/ageingM
#SBATCH --output=logs/influenza_ts_fva91_%j.log

set -euo pipefail

echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Working directory: $(pwd)"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"

export OPENBLAS_NUM_THREADS=1
export JULIA_CPU_TARGET=generic
export JULIA_NUM_PRECOMPILE_TASKS=1

echo "Julia version:"
"$HOME/.juliaup/bin/julia" --version

echo "Preparing project dependencies..."
srun "$HOME/.juliaup/bin/julia" \
    --startup-file=no \
    --project=. \
    --threads="${SLURM_CPUS_PER_TASK}" \
    -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "Starting time-series simulation at: $(date)"
srun "$HOME/.juliaup/bin/julia" \
    --startup-file=no \
    --project=. \
    --threads="${SLURM_CPUS_PER_TASK}" \
    scripts/runs/ageing_influenza_timeseries.jl

echo "Simulation finished at: $(date)"
echo "Expected CSV outputs:"
echo "  output/time_series_figure4_baseline_fva91_R0_1p5.csv"
echo "  output/time_age_figure4_baseline_fva91_R0_1p5.csv"
