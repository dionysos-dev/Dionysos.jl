#!/bin/bash
#SBATCH --job-name=biped4d_merge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/merge_%j.out
#SBATCH --error=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/merge_%j.err

# Gathers the chunks into one file, one chunk after the other: the memory is that of a chunk, not of the
# abstraction (a few GB are enough, unlike the solve job). The number of chunks is read from the names of
# the chunk files (DIONYSOS_NCHUNKS overrides it). Submit it from the root of the repository, with the same
# DIONYSOS_TRANSITION_OUTDIR and DIONYSOS_ABSTRACTION_FILE as the chunks if you changed them:
#
#   sbatch problems/BipedRobot/4D_model_underactuated/slurm_array/submit_merge_chunks.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
DIR=problems/BipedRobot/4D_model_underactuated/slurm_array
mkdir -p "$DIR/logs"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GKSwstype=100

# The same directory as the chunks (see submit_compute_chunks.sbatch), and the output file, e.g.
#   export DIONYSOS_ABSTRACTION_FILE="$GLOBALSCRATCH/biped4d/normal_form_abstraction.jld2"
export DIONYSOS_TRANSITION_OUTDIR="${DIONYSOS_TRANSITION_OUTDIR:-$PWD/$DIR/out/transitions}"
export DIONYSOS_ABSTRACTION_FILE="${DIONYSOS_ABSTRACTION_FILE:-$PWD/$DIR/out/normal_form_abstraction.jld2}"

echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "DIONYSOS_NCHUNKS=${DIONYSOS_NCHUNKS:-(from the chunk files)}"
echo "DIONYSOS_TRANSITION_OUTDIR=$DIONYSOS_TRANSITION_OUTDIR"
echo "DIONYSOS_ABSTRACTION_FILE=$DIONYSOS_ABSTRACTION_FILE"

set +e
srun --export=ALL julia --project=problems/BipedRobot "$DIR/merge_chunks.jl"
STATUS=$?
set -e
echo "[merge] julia exited with status $STATUS at $(date '+%Y-%m-%d %H:%M:%S')"
if [ "$STATUS" -ne 0 ]; then
    echo "[merge] FAILED (status $STATUS): see the .out and .err logs; 'sacct -j $SLURM_JOB_ID --format=JobID,State,ExitCode,MaxRSS,Elapsed,ReqMem' says OUT_OF_MEMORY or TIMEOUT if the job was killed" >&2
fi
exit $STATUS
