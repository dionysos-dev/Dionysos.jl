#!/bin/bash
# ============================================================================
# Local test runner (no SLURM needed).
# Runs all four parallelism modes sequentially for a given Dionysos example,
# collecting JSON metrics for each mode. Optionally generates a LaTeX report.
#
# Usage:  bash parallel_tests/run_local.sh <example.jl> [NPARTS] [--report]
#
# Examples:
#   bash parallel_tests/run_local.sh utils/example_distributed.jl
#   bash parallel_tests/run_local.sh utils/example_distributed.jl 4 --report
#   bash parallel_tests/run_local.sh problems/BipedRobot/6D_model/robot_example.jl 4 --report
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
    echo "Usage: $0 <example.jl> [NPARTS] [--report]"
    echo ""
    echo "  <example.jl>  Path to the Julia example script (relative to project root)"
    echo "  [NPARTS]      Number of distributed partitions (default: 4)"
    echo "  --report      Generate a LaTeX/PDF comparison report after all runs"
    exit 1
fi

EXAMPLE_SCRIPT="${POSITIONAL[0]}"
NPARTS="${POSITIONAL[1]:-4}"

# Resolve project root (parent of parallel_tests/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WRAPPER="${SCRIPT_DIR}/src/run_and_collect.jl"

# Derive a short name from the example script for output directories
EXAMPLE_NAME="$(basename "${EXAMPLE_SCRIPT}" .jl)"
RESULTS_BASE="${SCRIPT_DIR}/results/${EXAMPLE_NAME}"

# Resolve the full path to the example
EXAMPLE_FULL="${PROJECT_DIR}/${EXAMPLE_SCRIPT}"
if [ ! -f "${EXAMPLE_FULL}" ]; then
    echo "ERROR: Example script not found: ${EXAMPLE_FULL}"
    exit 1
fi

# Detect CPU count
NCPUS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo "============================================================"
echo "  Dionysos Parallel Test Runner (local)"
echo "  Example:  ${EXAMPLE_SCRIPT}"
echo "  NPARTS:   ${NPARTS}"
echo "  CPUs:     ${NCPUS}"
echo "  Report:   ${GENERATE_REPORT}"
echo "  Project:  ${PROJECT_DIR}"
echo "============================================================"
echo ""

# ── 1/4  Serial ──────────────────────────────────────────────
OUTDIR="${RESULTS_BASE}/serial"
mkdir -p "${OUTDIR}"
echo "──── 1/4  Serial ────"
echo "  Output: ${OUTDIR}"
DIONYSOS_DISTRIBUTED=false \
DIONYSOS_THREADED=false \
DIONYSOS_NPARTS=1 \
DIONYSOS_OUTDIR="${OUTDIR}" \
julia --project="${PROJECT_DIR}" \
      "${WRAPPER}" "${EXAMPLE_FULL}" 2>&1 | tee "${OUTDIR}/run.log"
echo ""

# ── 2/4  Threaded ────────────────────────────────────────────
OUTDIR="${RESULTS_BASE}/threaded"
mkdir -p "${OUTDIR}"
echo "──── 2/4  Threaded (${NCPUS} threads) ────"
echo "  Output: ${OUTDIR}"
DIONYSOS_DISTRIBUTED=false \
DIONYSOS_THREADED=true \
DIONYSOS_NPARTS=1 \
DIONYSOS_OUTDIR="${OUTDIR}" \
julia --project="${PROJECT_DIR}" \
      --threads="${NCPUS}" \
      "${WRAPPER}" "${EXAMPLE_FULL}" 2>&1 | tee "${OUTDIR}/run.log"
echo ""

# ── 3/4  Distributed ────────────────────────────────────────
OUTDIR="${RESULTS_BASE}/distributed"
mkdir -p "${OUTDIR}"
echo "──── 3/4  Distributed (${NPARTS} workers) ────"
echo "  Output: ${OUTDIR}"
DIONYSOS_DISTRIBUTED=true \
DIONYSOS_THREADED=false \
DIONYSOS_NPARTS="${NPARTS}" \
DIONYSOS_OUTDIR="${OUTDIR}" \
julia --project="${PROJECT_DIR}" \
      "${WRAPPER}" "${EXAMPLE_FULL}" 2>&1 | tee "${OUTDIR}/run.log"
echo ""

# ── 4/4  Hybrid ──────────────────────────────────────────────
OUTDIR="${RESULTS_BASE}/hybrid"
mkdir -p "${OUTDIR}"
echo "──── 4/4  Hybrid (${NPARTS} workers × ${NCPUS} threads) ────"
echo "  Output: ${OUTDIR}"
DIONYSOS_DISTRIBUTED=true \
DIONYSOS_THREADED=true \
DIONYSOS_NPARTS="${NPARTS}" \
DIONYSOS_OUTDIR="${OUTDIR}" \
julia --project="${PROJECT_DIR}" \
      --threads="${NCPUS}" \
      "${WRAPPER}" "${EXAMPLE_FULL}" 2>&1 | tee "${OUTDIR}/run.log"
echo ""

# ── Report generation ────────────────────────────────────────
if [ "${GENERATE_REPORT}" = "true" ]; then
    REPORT_DIR="${SCRIPT_DIR}/report/${EXAMPLE_NAME}"
    echo "──── Generating report ────"
    julia --project="${PROJECT_DIR}" \
          "${SCRIPT_DIR}/src/generate_report.jl" \
          "${RESULTS_BASE}" "${REPORT_DIR}"
    echo ""
fi

echo "============================================================"
echo "  All done!"
echo "  Results: ${RESULTS_BASE}/"
if [ "${GENERATE_REPORT}" = "true" ]; then
    echo "  Report:  ${SCRIPT_DIR}/report/${EXAMPLE_NAME}/"
fi
echo "============================================================"
