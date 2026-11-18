#!/bin/bash
#SBATCH --job-name=influenza_ts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --chdir=/hpc/dla_lti/jvanmeenen/ageingM
#SBATCH --output=logs/influenza_ts_%j.log

set -euo pipefail

# Usage: sbatch test_ageing_influenza_timeseries.sh <R0> <nu_per_year> <first_age_s> [output_tag]
# Example: sbatch test_ageing_influenza_timeseries.sh 1.5 0.5 91
# nu = 1 / mean vaccine immunity duration in years.
if [[ $# -ne 3 && $# -ne 4 ]]; then
    echo "Usage: sbatch test_ageing_influenza_timeseries.sh <R0> <nu_per_year> <first_age_s> [output_tag]" >&2
    exit 2
fi

echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Working directory: $(pwd)"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Parameters: R0=$1, nu=$2/yr, first vaccination age=$3; tag=${4:-figure4_baseline}"

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
    scripts/runs/ageing_influenza_timeseries.jl "$@"

echo "Simulation finished at: $(date)"
