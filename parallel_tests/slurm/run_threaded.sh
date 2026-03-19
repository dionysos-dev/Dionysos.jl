#!/bin/bash
#SBATCH --job-name=dionysos_threaded
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=32G

# ── Configuration ──────────────────────────────────────────────
# Usage: sbatch slurm/run_threaded.sh <example.jl>
EXAMPLE_SCRIPT="${1:?Usage: sbatch run_threaded.sh <example.jl>}"
NTHREADS=${SLURM_CPUS_PER_TASK:-8}

# PROJECT_DIR is either exported by run_all.sh or derived from SLURM_SUBMIT_DIR
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-.}}"
EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
OUTDIR="${PROJECT_DIR}/parallel_tests/results/${EXAMPLE_NAME}/threaded"

export DIONYSOS_DISTRIBUTED=false
export DIONYSOS_THREADED=true
export DIONYSOS_NPARTS=1
export DIONYSOS_OUTDIR="${OUTDIR}"
export JULIA_NUM_THREADS="${NTHREADS}"

mkdir -p "${OUTDIR}" "${PROJECT_DIR}/parallel_tests/logs"

echo "=========================================="
echo "  SLURM Threaded – Dionysos"
echo "  Job ID:    ${SLURM_JOB_ID}"
echo "  Node:      $(hostname)"
echo "  Threads:   ${NTHREADS}"
echo "  Example:   ${EXAMPLE_SCRIPT}"
echo "=========================================="

WRAPPER="${PROJECT_DIR}/parallel_tests/src/run_and_collect.jl"
srun julia --project="${PROJECT_DIR}" \
     --threads="${NTHREADS}" \
     "${WRAPPER}" "${PROJECT_DIR}/${EXAMPLE_SCRIPT}" 2>&1 | tee "${OUTDIR}/run.log"

echo "Job finished at $(date)"
