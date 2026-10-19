#!/bin/bash
#SBATCH --job-name=biped4d_precompile
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/precompile_%j.out
#SBATCH --error=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/precompile_%j.err

# Builds the Julia precompilation cache once, on a compute node, before the array (`submit_all.sh` submits
# it and makes the array wait for it): without it, every task of the array compiles the packages itself,
# for minutes, and the 1-hour budget of a chunk goes there. It must run on the same kind of node as the
# tasks: give it the same partition/constraints as `submit_compute_chunks.sbatch` (add the same
# `#SBATCH --partition=...` / `--constraint=...` lines to both). Submit it alone, from the root of the
# repository, to rebuild the cache: sbatch problems/BipedRobot/4D_model_underactuated/slurm_array/submit_precompile.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
DIR=problems/BipedRobot/4D_model_underactuated/slurm_array
mkdir -p "$DIR/logs"

export JULIA_NUM_THREADS=1
export JULIA_NUM_PRECOMPILE_TASKS="${SLURM_CPUS_PER_TASK:-1}"   # the cores of the job, not those of the node
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GKSwstype=100

echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "JULIA_CPU_TARGET=${JULIA_CPU_TARGET:-(default)}"

set +e
srun --export=ALL julia --project=problems/BipedRobot "$DIR/precompile.jl"
STATUS=$?
set -e
echo "[precompile] julia exited with status $STATUS at $(date '+%Y-%m-%d %H:%M:%S')"
if [ "$STATUS" -ne 0 ]; then
    echo "[precompile] FAILED (status $STATUS): see the .out and .err logs" >&2
fi
exit $STATUS
