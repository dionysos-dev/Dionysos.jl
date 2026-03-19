#!/bin/bash
#SBATCH --job-name=dionysos_bench
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=8G

# ============================================================================
# Generic SLURM runner for any Dionysos example in any parallelism mode.
#
# Usage:
#   sbatch [SLURM_OPTS] run_example.sh <example.jl> <distributed> <threaded>
#
#   <example.jl>    Path to Julia script relative to project root
#   <distributed>   true/false – enable distributed workers
#   <threaded>      true/false – enable multi-threading
#
# When using distributed mode, set --ntasks > 1 in your sbatch options.
# When using threaded mode, set --cpus-per-task > 1 in your sbatch options.
#
# Examples:
#   # Serial
#   sbatch --ntasks=1 --cpus-per-task=1 run_example.sh utils/example_distributed.jl false false
#
#   # Threaded
#   sbatch --ntasks=1 --cpus-per-task=8 run_example.sh utils/example_distributed.jl false true
#
#   # Distributed
#   sbatch --ntasks=8 --cpus-per-task=1 run_example.sh utils/example_distributed.jl true false
#
#   # Hybrid
#   sbatch --ntasks=4 --cpus-per-task=4 run_example.sh utils/example_distributed.jl true true
# ============================================================================

EXAMPLE_SCRIPT="${1:?Usage: sbatch run_example.sh <example.jl> <distributed> <threaded>}"
DISTRIBUTED="${2:?Specify distributed: true or false}"
THREADED="${3:?Specify threaded: true or false}"

# Derive mode tag
if [ "${DISTRIBUTED}" = "true" ] && [ "${THREADED}" = "true" ]; then
    MODE="hybrid"
elif [ "${DISTRIBUTED}" = "true" ]; then
    MODE="distributed"
elif [ "${THREADED}" = "true" ]; then
    MODE="threaded"
else
    MODE="serial"
fi

# Project directory
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-.}}"
EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
OUTDIR="${PROJECT_DIR}/parallel_tests/results/${EXAMPLE_NAME}/${MODE}"

# Threading
NTHREADS=${SLURM_CPUS_PER_TASK:-1}
NPARTS=${SLURM_NTASKS:-1}

export DIONYSOS_DISTRIBUTED="${DISTRIBUTED}"
export DIONYSOS_THREADED="${THREADED}"
export DIONYSOS_NPARTS="${NPARTS}"
export DIONYSOS_OUTDIR="${OUTDIR}"
export JULIA_NUM_THREADS="${NTHREADS}"

mkdir -p "${OUTDIR}" "${PROJECT_DIR}/parallel_tests/logs"

echo "=========================================="
echo "  Dionysos SLURM Runner"
echo "  Mode:      ${MODE}"
echo "  Job ID:    ${SLURM_JOB_ID:-local}"
echo "  Node:      $(hostname)"
echo "  Tasks:     ${NPARTS}"
echo "  Threads:   ${NTHREADS} per task"
echo "  Example:   ${EXAMPLE_SCRIPT}"
echo "  Output:    ${OUTDIR}"
echo "=========================================="

# Build julia command
WRAPPER="${PROJECT_DIR}/parallel_tests/src/run_and_collect.jl"
JULIA_CMD="julia --project=${PROJECT_DIR}"
if [ "${THREADED}" = "true" ] && [ "${NTHREADS}" -gt 1 ]; then
    JULIA_CMD="${JULIA_CMD} --threads=${NTHREADS}"
fi

if [ -n "${SLURM_JOB_ID:-}" ]; then
    srun ${JULIA_CMD} "${WRAPPER}" "${PROJECT_DIR}/${EXAMPLE_SCRIPT}" 2>&1 | tee "${OUTDIR}/run.log"
else
    ${JULIA_CMD} "${WRAPPER}" "${PROJECT_DIR}/${EXAMPLE_SCRIPT}" 2>&1 | tee "${OUTDIR}/run.log"
fi

echo "Job finished at $(date)"
