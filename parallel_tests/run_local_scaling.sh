#!/bin/bash
# ============================================================================
# Local scaling test runner.
# Sweeps a range of worker/thread counts for a given mode (distributed or
# threaded) while keeping other parameters fixed.
#
# Usage:
#   bash parallel_tests/run_local_scaling.sh <example.jl> <mode> <min> <max> [step] [--report]
#
# Arguments:
#   <example.jl>  Path to Julia script (relative to project root)
#   <mode>        One of: distributed, threaded
#   <min>         Minimum number of workers/threads
#   <max>         Maximum number of workers/threads
#   [step]        Step size (default: 1)
#   --report      Generate scaling report after all runs
#
# Examples:
#   bash parallel_tests/run_local_scaling.sh problems/BipedRobot/6D_model/robot_example.jl distributed 1 8
#   bash parallel_tests/run_local_scaling.sh problems/BipedRobot/6D_model/robot_example.jl threaded 1 8 2 --report
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

if [ ${#POSITIONAL[@]} -lt 4 ]; then
    echo "Usage: $0 <example.jl> <mode> <min> <max> [step] [--report]"
    echo ""
    echo "  <mode>   distributed | threaded"
    echo "  <min>    Minimum workers/threads"
    echo "  <max>    Maximum workers/threads"
    echo "  [step]   Step size (default: 1)"
    exit 1
fi

EXAMPLE_SCRIPT="${POSITIONAL[0]}"
MODE="${POSITIONAL[1]}"
N_MIN="${POSITIONAL[2]}"
N_MAX="${POSITIONAL[3]}"
N_STEP="${POSITIONAL[4]:-1}"

if [[ "$MODE" != "distributed" && "$MODE" != "threaded" ]]; then
    echo "ERROR: mode must be 'distributed' or 'threaded', got '$MODE'"
    exit 1
fi

# Resolve directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PARALLEL_ENV="${SCRIPT_DIR}"
WRAPPER="${SCRIPT_DIR}/src/run_and_collect.jl"

EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
RESULTS_BASE="${SCRIPT_DIR}/results/${EXAMPLE_NAME}"

EXAMPLE_FULL="${PROJECT_DIR}/${EXAMPLE_SCRIPT}"
if [ ! -f "${EXAMPLE_FULL}" ]; then
    echo "ERROR: Example script not found: ${EXAMPLE_FULL}"
    exit 1
fi

echo "============================================================"
echo "  Dionysos Scaling Test Runner (local)"
echo "  Example:  ${EXAMPLE_SCRIPT}"
echo "  Mode:     ${MODE}"
echo "  Range:    ${N_MIN} to ${N_MAX} (step ${N_STEP})"
echo "  Report:   ${GENERATE_REPORT}"
echo "  Project:  ${PROJECT_DIR}"
echo "============================================================"
echo ""

RUN_IDX=0
for N in $(seq "${N_MIN}" "${N_STEP}" "${N_MAX}"); do
    RUN_IDX=$((RUN_IDX + 1))
    OUTDIR="${RESULTS_BASE}/${MODE}_n${N}"
    mkdir -p "${OUTDIR}"

    if [ "$MODE" = "distributed" ]; then
        echo "──── Run ${RUN_IDX}: Distributed with ${N} workers ────"
        echo "  Output: ${OUTDIR}"
        DIONYSOS_DISTRIBUTED=true \
        DIONYSOS_THREADED=false \
        DIONYSOS_NPARTS="${N}" \
        DIONYSOS_OUTDIR="${OUTDIR}" \
        julia --project="${PARALLEL_ENV}" \
              "${WRAPPER}" "${EXAMPLE_FULL}" 2>&1 | tee "${OUTDIR}/run.log"
    else
        echo "──── Run ${RUN_IDX}: Threaded with ${N} threads ────"
        echo "  Output: ${OUTDIR}"
        DIONYSOS_DISTRIBUTED=false \
        DIONYSOS_THREADED=true \
        DIONYSOS_NPARTS=1 \
        DIONYSOS_OUTDIR="${OUTDIR}" \
        julia --project="${PARALLEL_ENV}" \
              --threads="${N}" \
              "${WRAPPER}" "${EXAMPLE_FULL}" 2>&1 | tee "${OUTDIR}/run.log"
    fi
    echo ""
done

# ── Report generation ────────────────────────────────────────
if [ "${GENERATE_REPORT}" = "true" ]; then
    REPORT_DIR="${SCRIPT_DIR}/report/${EXAMPLE_NAME}"
    echo "──── Generating report ────"
    julia --project="${PARALLEL_ENV}" \
          "${SCRIPT_DIR}/src/generate_report.jl" \
          "${RESULTS_BASE}" "${REPORT_DIR}"
    echo ""
fi

echo "============================================================"
echo "  Scaling sweep done!"
echo "  Results: ${RESULTS_BASE}/"
if [ "${GENERATE_REPORT}" = "true" ]; then
    echo "  Report:  ${SCRIPT_DIR}/report/${EXAMPLE_NAME}/"
fi
echo "============================================================"
