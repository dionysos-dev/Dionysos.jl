#!/bin/bash
#SBATCH --job-name=biped4d_solve
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=02:30:00
#SBATCH --output=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/solve_%j.out
#SBATCH --error=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/solve_%j.err

# Control synthesis, closed loop and animation on the merged abstraction. Submitted after the merge by
# `submit_all.sh`; by hand, from the root of the repository, to rerun it with another start or tolerance:
#
#   DIONYSOS_SOLVE_PARAMS=my_params.jl DIONYSOS_SCENARIO=tol3cm \
#       sbatch problems/BipedRobot/4D_model_underactuated/slurm_array/submit_solve.sbatch
#
# The memory is that of the merge job (the abstraction is loaded into the automaton again) plus the
# synthesis: set it from the merge job (`seff <jobid>`), with a margin. Everything else is in solve.jl.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
DIR=problems/BipedRobot/4D_model_underactuated/slurm_array
mkdir -p "$DIR/logs"
source "$DIR/env.sh"     # the precompilation cache of this CPU type

export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GKSwstype=100                      # headless Plots: the animation is written to files

# The merged file (the same variable as the merge job) and where the figures go, e.g.
#   export DIONYSOS_ABSTRACTION_FILE="$GLOBALSCRATCH/biped4d/normal_form_abstraction.jld2"
export DIONYSOS_ABSTRACTION_FILE="${DIONYSOS_ABSTRACTION_FILE:-$PWD/$DIR/out/normal_form_abstraction.jld2}"
export DIONYSOS_OUTPUT_DIR="${DIONYSOS_OUTPUT_DIR:-$PWD/$DIR/out/results}"
export DIONYSOS_SCENARIO="${DIONYSOS_SCENARIO:-job${SLURM_JOB_ID}}"
if [ -n "${DIONYSOS_SOLVE_PARAMS:-}" ]; then DIONYSOS_SOLVE_PARAMS="$(realpath "$DIONYSOS_SOLVE_PARAMS")"; export DIONYSOS_SOLVE_PARAMS; fi

echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "DIONYSOS_ABSTRACTION_FILE=$DIONYSOS_ABSTRACTION_FILE"
echo "DIONYSOS_OUTPUT_DIR=$DIONYSOS_OUTPUT_DIR"
echo "DIONYSOS_SCENARIO=$DIONYSOS_SCENARIO"
echo "DIONYSOS_SOLVE_PARAMS=${DIONYSOS_SOLVE_PARAMS:-}"

set +e
srun --export=ALL julia --project=problems/BipedRobot \
    --heap-size-hint="$(( ${SLURM_MEM_PER_NODE:-65536} * 9 / 10 ))M" \
    "$DIR/solve.jl"
STATUS=$?
set -e
echo "[solve] julia exited with status $STATUS at $(date '+%Y-%m-%d %H:%M:%S')"
if [ "$STATUS" -ne 0 ]; then
    echo "[solve] FAILED (status $STATUS): see the .out and .err logs; 'sacct -j $SLURM_JOB_ID --format=JobID,State,ExitCode,MaxRSS,Elapsed,ReqMem' says OUT_OF_MEMORY or TIMEOUT if the job was killed" >&2
fi
exit $STATUS
