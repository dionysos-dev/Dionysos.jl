#!/bin/bash
#SBATCH --job-name=dionysos_scale_thr
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=32G

# ============================================================================
# SLURM array job for threaded scaling.
# Each array task runs with a different number of threads (SLURM_ARRAY_TASK_ID).
#
# Usage:
#   sbatch --array=1-8 slurm/run_scaling_threaded.sh <example.jl>
#   sbatch --array=1,2,4,8 slurm/run_scaling_threaded.sh <example.jl>
#
# The array index is used as the number of threads.
# --cpus-per-task is set dynamically via the wrapper below.
# ============================================================================
set -euo pipefail

EXAMPLE_SCRIPT="${1:?Usage: sbatch --array=1-8 run_scaling_threaded.sh <example.jl>}"
NTHREADS="${SLURM_ARRAY_TASK_ID}"

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-.}}"
PARALLEL_ENV="${PROJECT_DIR}/parallel_tests"
EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
OUTDIR="${PROJECT_DIR}/parallel_tests/results/${EXAMPLE_NAME}/threaded_n${NTHREADS}"

export DIONYSOS_DISTRIBUTED=false
export DIONYSOS_THREADED=true
export DIONYSOS_NPARTS=1
export DIONYSOS_OUTDIR="${OUTDIR}"
export JULIA_NUM_THREADS="${NTHREADS}"

mkdir -p "${OUTDIR}" "${PROJECT_DIR}/parallel_tests/logs"

echo "=========================================="
echo "  SLURM Threaded Scaling – Dionysos"
echo "  Job ID:    ${SLURM_JOB_ID} (array ${SLURM_ARRAY_TASK_ID})"
echo "  Threads:   ${NTHREADS}"
echo "  Node:      $(hostname)"
echo "  Example:   ${EXAMPLE_SCRIPT}"
echo "=========================================="

WRAPPER="${PROJECT_DIR}/parallel_tests/src/run_and_collect.jl"
srun julia --project="${PARALLEL_ENV}" \
     --threads="${NTHREADS}" \
     "${WRAPPER}" "${PROJECT_DIR}/${EXAMPLE_SCRIPT}" 2>&1 | tee "${OUTDIR}/run.log"

echo "Job finished at $(date)"
