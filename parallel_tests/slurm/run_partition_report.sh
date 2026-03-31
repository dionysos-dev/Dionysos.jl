#!/bin/bash
# ==============================================================================
#  SLURM job — generate a LaTeX report for a partition pipeline run.
#
#  Environment variables:
#      DIONYSOS_PROJECT_ROOT — absolute path to the Dionysos repository root
#      DIONYSOS_PARTDIR      — directory containing partition_*.jld2 files
#      DIONYSOS_REPORTDIR    — output directory for report files (default: PARTDIR/report)
# ==============================================================================

#SBATCH --job-name=dionysos_part_report
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --output=logs/partition_report_%j.out
#SBATCH --error=logs/partition_report_%j.err

set -euo pipefail

if [ -n "${DIONYSOS_PROJECT_ROOT:-}" ]; then
    PROJECT_ROOT="${DIONYSOS_PROJECT_ROOT}"
elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "${SLURM_SUBMIT_DIR}/parallel_tests" ]; then
    PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi

PARTITIONS_DIR="${DIONYSOS_PARTDIR:?DIONYSOS_PARTDIR must be set}"
REPORT_DIR="${DIONYSOS_REPORTDIR:-${PARTITIONS_DIR}/report}"
PARALLEL_ENV="${PROJECT_ROOT}/parallel_tests"

echo "============================================================"
echo "  SLURM Partition Report Job"
echo "============================================================"
echo "  Job ID        : ${SLURM_JOB_ID}"
echo "  Project root  : ${PROJECT_ROOT}"
echo "  Partitions dir: ${PARTITIONS_DIR}"
echo "  Report dir    : ${REPORT_DIR}"
echo "============================================================"
echo

srun julia --project="${PARALLEL_ENV}" \
    "${PROJECT_ROOT}/parallel_tests/src/generate_partition_report.jl" \
    "${PARTITIONS_DIR}" \
    "${REPORT_DIR}"
