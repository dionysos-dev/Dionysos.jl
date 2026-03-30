#!/bin/bash
# ==============================================================================
#  Submit the full SLURM array pipeline:
#    1. Array job  — one task per partition (fast allocation, ~5-15 min each)
#    2. Assembly job — merges all partition results (runs after array completes)
#
#  Usage:
#      bash slurm/submit_array_pipeline.sh [OPTIONS]
#
#  Options:
#      -n, --nparts NUM        Number of partitions (default: 300)
#      -s, --strategy STR      Partition strategy: roundrobin|contiguous (default: roundrobin)
#      -t, --time TIME         Wall-time per partition task (default: 00:15:00)
#      -m, --mem MEM           Memory per partition task (default: 8G)
#      --assemble-time TIME    Wall-time for assembly job (default: 00:30:00)
#      --assemble-mem MEM      Memory for assembly job (default: 16G)
#      --solve                 Also solve optimal-control after assembly
#      --setup PATH            Path to setup script
#      --outdir DIR            Output directory (default: auto)
#      --dry-run               Print commands but do not submit
# ==============================================================================

set -euo pipefail

# --- Defaults ---
NPARTS=300
STRATEGY="roundrobin"
PART_TIME="00:15:00"
PART_MEM="8G"
ASSEMBLE_TIME="00:30:00"
ASSEMBLE_MEM="16G"
SOLVE="false"
DRY_RUN=false

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SETUP_SCRIPT="${PROJECT_ROOT}/problems/BipedRobot/6D_model/robot_example_setup.jl"
OUTDIR=""

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--nparts)        NPARTS="$2";        shift 2 ;;
        -s|--strategy)      STRATEGY="$2";      shift 2 ;;
        -t|--time)          PART_TIME="$2";     shift 2 ;;
        -m|--mem)           PART_MEM="$2";      shift 2 ;;
        --assemble-time)    ASSEMBLE_TIME="$2"; shift 2 ;;
        --assemble-mem)     ASSEMBLE_MEM="$2";  shift 2 ;;
        --solve)            SOLVE="true";       shift   ;;
        --setup)            SETUP_SCRIPT="$2";  shift 2 ;;
        --outdir)           OUTDIR="$2";        shift 2 ;;
        --dry-run)          DRY_RUN=true;       shift   ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# --- Determine output directory ---
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
if [ -z "${OUTDIR}" ]; then
    OUTDIR="${PROJECT_ROOT}/parallel_tests/results/robot_example/array_${TIMESTAMP}"
fi

# --- Create directories ---
mkdir -p "${PROJECT_ROOT}/parallel_tests/slurm/logs"
mkdir -p "${OUTDIR}"

echo "============================================================"
echo "  Dionysos SLURM Array Pipeline"
echo "============================================================"
echo "  N parts       : ${NPARTS}"
echo "  Strategy      : ${STRATEGY}"
echo "  Part time     : ${PART_TIME}"
echo "  Part memory   : ${PART_MEM}"
echo "  Assemble time : ${ASSEMBLE_TIME}"
echo "  Assemble mem  : ${ASSEMBLE_MEM}"
echo "  Solve         : ${SOLVE}"
echo "  Setup script  : ${SETUP_SCRIPT}"
echo "  Output dir    : ${OUTDIR}"
echo "  Dry run       : ${DRY_RUN}"
echo "============================================================"
echo

# --- Step 1: Submit the array job ---
ARRAY_CMD="sbatch --parsable \
    --array=1-${NPARTS} \
    --time=${PART_TIME} \
    --mem=${PART_MEM} \
    --output=${PROJECT_ROOT}/parallel_tests/slurm/logs/partition_%A_%a.out \
    --error=${PROJECT_ROOT}/parallel_tests/slurm/logs/partition_%A_%a.err \
    --export=ALL,DIONYSOS_PROJECT_ROOT=${PROJECT_ROOT},DIONYSOS_NPARTS=${NPARTS},DIONYSOS_STRATEGY=${STRATEGY},DIONYSOS_SETUP=${SETUP_SCRIPT},DIONYSOS_OUTDIR=${OUTDIR} \
    ${SCRIPT_DIR}/run_array_partitions.sh"

echo "Step 1: Submitting array job (${NPARTS} tasks)..."
echo "  ${ARRAY_CMD}"
echo

if [ "${DRY_RUN}" = true ]; then
    ARRAY_JOB_ID="DRYRUN_123456"
    echo "  [DRY RUN] Would submit array job -> ${ARRAY_JOB_ID}"
else
    ARRAY_JOB_ID=$(eval ${ARRAY_CMD})
    echo "  Array job submitted: ${ARRAY_JOB_ID}"
fi

# --- Step 2: Submit the assembly job with dependency ---
ASSEMBLE_CMD="sbatch --parsable \
    --dependency=afterok:${ARRAY_JOB_ID} \
    --time=${ASSEMBLE_TIME} \
    --mem=${ASSEMBLE_MEM} \
    --output=${PROJECT_ROOT}/parallel_tests/slurm/logs/assemble_%j.out \
    --error=${PROJECT_ROOT}/parallel_tests/slurm/logs/assemble_%j.err \
    --export=ALL,DIONYSOS_PROJECT_ROOT=${PROJECT_ROOT},DIONYSOS_NPARTS=${NPARTS},DIONYSOS_STRATEGY=${STRATEGY},DIONYSOS_SETUP=${SETUP_SCRIPT},DIONYSOS_PARTDIR=${OUTDIR},DIONYSOS_OUTDIR=${OUTDIR},DIONYSOS_SOLVE=${SOLVE} \
    ${SCRIPT_DIR}/run_assemble.sh"

echo
echo "Step 2: Submitting assembly job (dependency: afterok:${ARRAY_JOB_ID})..."
echo "  ${ASSEMBLE_CMD}"
echo

if [ "${DRY_RUN}" = true ]; then
    ASSEMBLE_JOB_ID="DRYRUN_789012"
    echo "  [DRY RUN] Would submit assembly job -> ${ASSEMBLE_JOB_ID}"
else
    ASSEMBLE_JOB_ID=$(eval ${ASSEMBLE_CMD})
    echo "  Assembly job submitted: ${ASSEMBLE_JOB_ID}"
fi

# --- Summary ---
echo
echo "============================================================"
echo "  Pipeline submitted!"
echo "============================================================"
echo "  Array job ID   : ${ARRAY_JOB_ID}  (${NPARTS} tasks)"
echo "  Assembly job ID: ${ASSEMBLE_JOB_ID}  (afterok:${ARRAY_JOB_ID})"
echo "  Output dir     : ${OUTDIR}"
echo
echo "  Monitor with:"
echo "    squeue -u \$USER"
echo "    sacct -j ${ARRAY_JOB_ID} --format=JobID,State,Elapsed,MaxRSS"
echo "    sacct -j ${ASSEMBLE_JOB_ID} --format=JobID,State,Elapsed,MaxRSS"
echo
echo "  Check partition logs:"
echo "    tail -f ${PROJECT_ROOT}/parallel_tests/slurm/logs/partition_${ARRAY_JOB_ID}_*.out"
echo
echo "  Check assembly log:"
echo "    tail -f ${PROJECT_ROOT}/parallel_tests/slurm/logs/assemble_${ASSEMBLE_JOB_ID}.out"
echo "============================================================"
