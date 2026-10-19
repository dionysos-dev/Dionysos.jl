#!/bin/bash
#SBATCH --job-name=biped4d_chunks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/chunk_%A_%a.out
#SBATCH --error=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/chunk_%A_%a.err

# One task of the array = one chunk of the abstraction, on one core. The array size is NOT set here:
# use `submit_all.sh NCHUNKS`, or `DIONYSOS_NCHUNKS=N sbatch --array=1-N <this file>` from the root of the
# repository (to re-run failed tasks only: `DIONYSOS_NCHUNKS=N sbatch --array=17,42 <this file>`).

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
DIR=problems/BipedRobot/4D_model_underactuated/slurm_array
mkdir -p "$DIR/logs"
source "$DIR/env.sh"     # the precompilation cache of this CPU type

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GKSwstype=100

# All the tasks write in the same directory; use a scratch directory for a big abstraction, e.g.
#   export DIONYSOS_TRANSITION_OUTDIR="$GLOBALSCRATCH/biped4d/transitions"
export DIONYSOS_TRANSITION_OUTDIR="${DIONYSOS_TRANSITION_OUTDIR:-$PWD/$DIR/out/transitions}"
mkdir -p "$DIONYSOS_TRANSITION_OUTDIR"

echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT"
echo "DIONYSOS_TRANSITION_OUTDIR=$DIONYSOS_TRANSITION_OUTDIR"

set +e
srun --export=ALL julia --project=problems/BipedRobot "$DIR/compute_chunk.jl"
STATUS=$?
set -e
echo "[chunk $SLURM_ARRAY_TASK_ID] julia exited with status $STATUS at $(date '+%Y-%m-%d %H:%M:%S')"
if [ "$STATUS" -ne 0 ]; then
    echo "[chunk $SLURM_ARRAY_TASK_ID] FAILED (status $STATUS): see the .out and .err logs; 'sacct -j $SLURM_JOB_ID --format=JobID,State,ExitCode,MaxRSS,Elapsed,ReqMem' says OUT_OF_MEMORY or TIMEOUT if the job was killed" >&2
fi
exit $STATUS
