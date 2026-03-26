#!/bin/bash
#SBATCH --job-name=dionysos_scale_dist
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=8G

# ============================================================================
# SLURM array job for distributed scaling.
# Each array task runs with a different number of workers (SLURM_ARRAY_TASK_ID).
#
# Usage:
#   sbatch --array=1-8 slurm/run_scaling_distributed.sh <example.jl>
#   sbatch --array=1,2,4,8 slurm/run_scaling_distributed.sh <example.jl>
#
# The array index is used as the number of workers (NPARTS).
# Adjust --ntasks-per-node or rely on Julia addprocs.
# ============================================================================
set -euo pipefail

EXAMPLE_SCRIPT="${1:?Usage: sbatch --array=1-8 run_scaling_distributed.sh <example.jl>}"
NWORKERS="${SLURM_ARRAY_TASK_ID}"

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-.}}"
PARALLEL_ENV="${PROJECT_DIR}/parallel_tests"
EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
OUTDIR="${PROJECT_DIR}/parallel_tests/results/${EXAMPLE_NAME}/distributed_n${NWORKERS}"

export DIONYSOS_DISTRIBUTED=true
export DIONYSOS_THREADED=false
export DIONYSOS_NPARTS="${NWORKERS}"
export DIONYSOS_OUTDIR="${OUTDIR}"

mkdir -p "${OUTDIR}" "${PROJECT_DIR}/parallel_tests/logs"

echo "=========================================="
echo "  SLURM Distributed Scaling – Dionysos"
echo "  Job ID:    ${SLURM_JOB_ID} (array ${SLURM_ARRAY_TASK_ID})"
echo "  Workers:   ${NWORKERS}"
echo "  Node:      $(hostname)"
echo "  Example:   ${EXAMPLE_SCRIPT}"
echo "=========================================="

WRAPPER="${PROJECT_DIR}/parallel_tests/src/run_and_collect.jl"
srun julia --project="${PARALLEL_ENV}" \
     "${WRAPPER}" "${PROJECT_DIR}/${EXAMPLE_SCRIPT}" 2>&1 | tee "${OUTDIR}/run.log"

echo "Job finished at $(date)"
