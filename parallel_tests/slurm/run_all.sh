#!/bin/bash
# ============================================================================
# Master script: submit all 4 parallelism modes to SLURM for a given example.
# Usage:  bash parallel_tests/slurm/run_all.sh <example.jl> [--report]
#
# Examples:
#   bash parallel_tests/slurm/run_all.sh utils/example_distributed.jl
#   bash parallel_tests/slurm/run_all.sh utils/example_distributed.jl --report
#   bash parallel_tests/slurm/run_all.sh problems/BipedRobot/6D_model/robot_example.jl --report
# ============================================================================
set -euo pipefail

# ── Parse arguments ──────────────────────────────────────────
GENERATE_REPORT=false
POSITIONAL=()
for arg in "$@"; do
    case "$arg" in
        --report) GENERATE_REPORT=true ;;
        *)        POSITIONAL+=("$arg") ;;
    esac
done

if [ ${#POSITIONAL[@]} -lt 1 ]; then
    echo "Usage: $0 <example.jl> [--report]"
    echo ""
    echo "  <example.jl>  Path to the Julia example script (relative to project root)"
    echo "  --report      Generate LaTeX/PDF report after all benchmarks finish"
    exit 1
fi

EXAMPLE_SCRIPT="${POSITIONAL[0]}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
RESULTS_DIR="${PROJECT_DIR}/parallel_tests/results/${EXAMPLE_NAME}"
REPORT_DIR="${PROJECT_DIR}/parallel_tests/report/${EXAMPLE_NAME}"

mkdir -p "${PROJECT_DIR}/parallel_tests/logs" "${RESULTS_DIR}"

echo "============================================================"
echo "  Dionysos Parallel Test Suite (SLURM)"
echo "  Example:  ${EXAMPLE_SCRIPT}"
echo "  Report:   ${GENERATE_REPORT}"
echo "  Project:  ${PROJECT_DIR}"
echo "============================================================"
echo ""

# Submit jobs and capture job IDs
echo "Submitting serial benchmark …"
JID_SERIAL=$(sbatch --parsable \
    --export=ALL,PROJECT_DIR="${PROJECT_DIR}" \
    --output="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_serial_%j.out" \
    --error="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_serial_%j.err" \
    "${SCRIPT_DIR}/run_serial.sh" "${EXAMPLE_SCRIPT}")
echo "  Job ID: ${JID_SERIAL}"

echo "Submitting threaded benchmark …"
JID_THREADED=$(sbatch --parsable \
    --export=ALL,PROJECT_DIR="${PROJECT_DIR}" \
    --output="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_threaded_%j.out" \
    --error="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_threaded_%j.err" \
    "${SCRIPT_DIR}/run_threaded.sh" "${EXAMPLE_SCRIPT}")
echo "  Job ID: ${JID_THREADED}"

echo "Submitting distributed benchmark …"
JID_DISTRIBUTED=$(sbatch --parsable \
    --export=ALL,PROJECT_DIR="${PROJECT_DIR}" \
    --output="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_distributed_%j.out" \
    --error="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_distributed_%j.err" \
    "${SCRIPT_DIR}/run_distributed.sh" "${EXAMPLE_SCRIPT}")
echo "  Job ID: ${JID_DISTRIBUTED}"

echo "Submitting hybrid benchmark …"
JID_HYBRID=$(sbatch --parsable \
    --export=ALL,PROJECT_DIR="${PROJECT_DIR}" \
    --output="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_hybrid_%j.out" \
    --error="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_hybrid_%j.err" \
    "${SCRIPT_DIR}/run_hybrid.sh" "${EXAMPLE_SCRIPT}")
echo "  Job ID: ${JID_HYBRID}"

echo ""
echo "All benchmarks submitted."
echo "  Serial:      ${JID_SERIAL}"
echo "  Threaded:    ${JID_THREADED}"
echo "  Distributed: ${JID_DISTRIBUTED}"
echo "  Hybrid:      ${JID_HYBRID}"

# ── Report generation (dependent job) ────────────────────────
if [ "${GENERATE_REPORT}" = "true" ]; then
    DEPS="${JID_SERIAL}:${JID_THREADED}:${JID_DISTRIBUTED}:${JID_HYBRID}"
    echo ""
    echo "Submitting report generation (after all benchmarks) …"
    JID_REPORT=$(sbatch --parsable \
        --dependency="afterany:${DEPS}" \
        --export=ALL,PROJECT_DIR="${PROJECT_DIR}" \
        --job-name="dionysos_report" \
        --output="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_report_%j.out" \
        --error="${PROJECT_DIR}/parallel_tests/logs/${EXAMPLE_NAME}_report_%j.err" \
        --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=4G --time=00:30:00 \
        --wrap="julia --project=\"${PROJECT_DIR}\" \
                \"${PROJECT_DIR}/parallel_tests/src/generate_report.jl\" \
                \"${RESULTS_DIR}\" \"${REPORT_DIR}\"")
    echo "  Report Job ID: ${JID_REPORT} (depends on ${DEPS})"
fi

echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Results will be in: ${RESULTS_DIR}/"
echo "Logs will be in:    ${PROJECT_DIR}/parallel_tests/logs/"
if [ "${GENERATE_REPORT}" = "true" ]; then
    echo "Report will be in:  ${REPORT_DIR}/"
fi
