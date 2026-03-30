#!/bin/bash
# ==============================================================================
#  SLURM job — assemble partition results into a complete abstraction.
#
#  Typically submitted with a dependency on the array job:
#      sbatch --dependency=afterok:<ARRAY_JOB_ID> slurm/run_assemble.sh
#
#  Or via the submit_array_pipeline.sh wrapper.
#
#  Environment variables (set before submission or via --export):
#      DIONYSOS_NPARTS      — total number of partitions
#      DIONYSOS_STRATEGY    — partition strategy (default: roundrobin)
#      DIONYSOS_SETUP       — path to setup script (default: auto-detected)
#      DIONYSOS_PARTDIR     — directory containing partition_*.jld2 files
#      DIONYSOS_OUTDIR      — output directory for assembled results (default: PARTDIR)
#      DIONYSOS_SOLVE       — if "true", also solve the optimal-control problem
# ==============================================================================

#SBATCH --job-name=dionysos_assemble
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --output=logs/assemble_%j.out
#SBATCH --error=logs/assemble_%j.err

set -euo pipefail

# --- Resolve paths ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PARALLEL_ENV="${PROJECT_ROOT}/parallel_tests"

NPARTS="${DIONYSOS_NPARTS:?DIONYSOS_NPARTS must be set}"
STRATEGY="${DIONYSOS_STRATEGY:-roundrobin}"
SETUP_SCRIPT="${DIONYSOS_SETUP:-${PROJECT_ROOT}/problems/BipedRobot/6D_model/robot_example_setup.jl}"
PARTITIONS_DIR="${DIONYSOS_PARTDIR:?DIONYSOS_PARTDIR must be set}"
OUTPUT_DIR="${DIONYSOS_OUTDIR:-${PARTITIONS_DIR}}"
SOLVE_FLAG=""
if [ "${DIONYSOS_SOLVE:-false}" = "true" ]; then
    SOLVE_FLAG="--solve"
fi

mkdir -p "$(dirname "${SLURM_OUTPUT_FILE:-logs/dummy}")" 2>/dev/null || true
mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "  SLURM Assembly Job"
echo "============================================================"
echo "  Job ID        : ${SLURM_JOB_ID}"
echo "  N parts       : ${NPARTS}"
echo "  Strategy      : ${STRATEGY}"
echo "  Node          : $(hostname)"
echo "  Partitions dir: ${PARTITIONS_DIR}"
echo "  Output dir    : ${OUTPUT_DIR}"
echo "  Solve         : ${DIONYSOS_SOLVE:-false}"
echo "============================================================"
echo

# --- Run the assembler ---
srun julia --project="${PARALLEL_ENV}" \
    "${PROJECT_ROOT}/parallel_tests/src/assemble_partitions.jl" \
    "${SETUP_SCRIPT}" \
    "${PARTITIONS_DIR}" \
    "${NPARTS}" \
    "${OUTPUT_DIR}" \
    "--strategy=${STRATEGY}" \
    ${SOLVE_FLAG}

echo
echo "Assembly finished at $(date)"
