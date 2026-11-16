#!/bin/bash
# Submit a Slurm array plus an automatic merge after all forcing tasks succeed.
# Usage: bash submit_influenza_array.sh 80  (or 30)
set -euo pipefail

coverage="${1:-}"
case "$coverage" in
    30|80) ;;
    *) echo "Usage: bash submit_influenza_array.sh <30|80>" >&2; exit 2 ;;
esac

cd "$(dirname "$0")"
mkdir -p logs output/forcing_parts
max_parallel="${INFLUENZA_MAX_PARALLEL:-3}"
[[ "$max_parallel" =~ ^[1-6]$ ]] || { echo "INFLUENZA_MAX_PARALLEL must be an integer from 1 to 6" >&2; exit 2; }

array_submission=$(sbatch --parsable \
    --job-name="influenza_c${coverage}_array" \
    --array="1-6%${max_parallel}" \
    --export="ALL,INFLUENZA_COVERAGE=${coverage}" \
    test_ageing_influenza_array.sh)
array_id="${array_submission%%;*}"

echo "Array job ID: ${array_id} (coverage ${coverage}%, max concurrent tasks: ${max_parallel})"
merge_submission=$(sbatch --parsable \
    --job-name="influenza_c${coverage}_merge" \
    --dependency="afterok:${array_id}" \
    --export="ALL,INFLUENZA_COVERAGE=${coverage},INFLUENZA_BATCH_ID=${array_id}" \
    test_merge_influenza_sweeps.sh)
merge_id="${merge_submission%%;*}"
echo "Merge job ID: ${merge_id} (runs after all array tasks succeed)"
echo "Results will be in: output/*_coverage${coverage}_batch${array_id}.csv"
echo "Progress: squeue -u \"\$USER\""
