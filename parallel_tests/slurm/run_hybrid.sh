#!/bin/bash
#SBATCH --job-name=dionysos_hybrid
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=8G

# ── Configuration ──────────────────────────────────────────────
# Usage: sbatch slurm/run_hybrid.sh <example.jl>
EXAMPLE_SCRIPT="${1:?Usage: sbatch run_hybrid.sh <example.jl>}"
NTHREADS=${SLURM_CPUS_PER_TASK:-4}
NPARTS=${SLURM_NTASKS:-4}

# PROJECT_DIR is either exported by run_all.sh or derived from SLURM_SUBMIT_DIR
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-.}}"
PARALLEL_ENV="${PROJECT_DIR}/parallel_tests"
EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
OUTDIR="${PROJECT_DIR}/parallel_tests/results/${EXAMPLE_NAME}/hybrid"

export DIONYSOS_DISTRIBUTED=true
export DIONYSOS_THREADED=true
export DIONYSOS_NPARTS="${NPARTS}"
export DIONYSOS_OUTDIR="${OUTDIR}"
export JULIA_NUM_THREADS="${NTHREADS}"

mkdir -p "${OUTDIR}" "${PROJECT_DIR}/parallel_tests/logs"

echo "=========================================="
echo "  SLURM Hybrid – Dionysos"
echo "  Job ID:    ${SLURM_JOB_ID}"
echo "  Nodes:     ${SLURM_NNODES}"
echo "  Tasks:     ${SLURM_NTASKS}"
echo "  Threads:   ${NTHREADS} per task"
echo "  Example:   ${EXAMPLE_SCRIPT}"
echo "=========================================="

WRAPPER="${PROJECT_DIR}/parallel_tests/src/run_and_collect.jl"
srun julia --project="${PARALLEL_ENV}" \
     --threads="${NTHREADS}" \
     "${WRAPPER}" "${PROJECT_DIR}/${EXAMPLE_SCRIPT}" 2>&1 | tee "${OUTDIR}/run.log"

echo "Job finished at $(date)"
