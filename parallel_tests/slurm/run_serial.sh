#!/bin/bash
#SBATCH --job-name=dionysos_serial
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G

# ── Configuration ──────────────────────────────────────────────
# Usage: sbatch slurm/run_serial.sh <example.jl>
#   example.jl is relative to the Dionysos project root.
EXAMPLE_SCRIPT="${1:?Usage: sbatch run_serial.sh <example.jl>}"

# PROJECT_DIR is either exported by run_all.sh or derived from SLURM_SUBMIT_DIR
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-.}}"
EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
OUTDIR="${PROJECT_DIR}/parallel_tests/results/${EXAMPLE_NAME}/serial"

export DIONYSOS_DISTRIBUTED=false
export DIONYSOS_THREADED=false
export DIONYSOS_NPARTS=1
export DIONYSOS_OUTDIR="${OUTDIR}"

mkdir -p "${OUTDIR}" "${PROJECT_DIR}/parallel_tests/logs"

echo "=========================================="
echo "  SLURM Serial – Dionysos"
echo "  Job ID:    ${SLURM_JOB_ID}"
echo "  Node:      $(hostname)"
echo "  Example:   ${EXAMPLE_SCRIPT}"
echo "=========================================="

WRAPPER="${PROJECT_DIR}/parallel_tests/src/run_and_collect.jl"
srun julia --project="${PROJECT_DIR}" \
     "${WRAPPER}" "${PROJECT_DIR}/${EXAMPLE_SCRIPT}" 2>&1 | tee "${OUTDIR}/run.log"

echo "Job finished at $(date)"
