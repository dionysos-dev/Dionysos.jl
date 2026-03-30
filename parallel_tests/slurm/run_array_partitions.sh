#!/bin/bash
# ==============================================================================
#  SLURM array job — compute transitions for each partition independently.
#
#  Submit with:
#      sbatch --array=1-300 slurm/run_array_partitions.sh
#
#  Or via the submit_array_pipeline.sh wrapper.
#
#  Environment variables (set before submission or via --export):
#      DIONYSOS_NPARTS      — total number of partitions (default: from array size)
#      DIONYSOS_STRATEGY    — partition strategy: roundrobin | contiguous (default: roundrobin)
#      DIONYSOS_SETUP       — path to setup script (default: auto-detected)
#      DIONYSOS_OUTDIR      — output directory for partition files (default: auto)
# ==============================================================================

#SBATCH --job-name=dionysos_part
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --output=logs/partition_%A_%a.out
#SBATCH --error=logs/partition_%A_%a.err

set -euo pipefail

# --- Resolve paths ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PARALLEL_ENV="${PROJECT_ROOT}/parallel_tests"

PARTITION_IDX="${SLURM_ARRAY_TASK_ID}"
NPARTS="${DIONYSOS_NPARTS:-${SLURM_ARRAY_TASK_COUNT:-${SLURM_ARRAY_TASK_MAX}}}"
STRATEGY="${DIONYSOS_STRATEGY:-roundrobin}"

SETUP_SCRIPT="${DIONYSOS_SETUP:-${PROJECT_ROOT}/problems/BipedRobot/6D_model/robot_example_setup.jl}"
OUTPUT_DIR="${DIONYSOS_OUTDIR:-${PROJECT_ROOT}/parallel_tests/results/robot_example/array_${SLURM_ARRAY_JOB_ID}}"

# Create log and output directories
mkdir -p "$(dirname "${SLURM_OUTPUT_FILE:-logs/dummy}")" 2>/dev/null || true
mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "  SLURM Array Partition Job"
echo "============================================================"
echo "  Job ID       : ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Partition    : ${PARTITION_IDX} / ${NPARTS}"
echo "  Strategy     : ${STRATEGY}"
echo "  Node         : $(hostname)"
echo "  Setup script : ${SETUP_SCRIPT}"
echo "  Output dir   : ${OUTPUT_DIR}"
echo "============================================================"
echo

# --- Run the partition runner ---
srun julia --project="${PARALLEL_ENV}" \
    "${PROJECT_ROOT}/parallel_tests/src/partition_runner.jl" \
    "${SETUP_SCRIPT}" \
    "${PARTITION_IDX}" \
    "${NPARTS}" \
    "${OUTPUT_DIR}" \
    "--strategy=${STRATEGY}"

echo
echo "Partition ${PARTITION_IDX} finished at $(date)"
