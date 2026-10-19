#!/bin/bash
#SBATCH --job-name=biped4d_merge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=250G
#SBATCH --time=02:00:00
#SBATCH --output=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/merge_%j.out
#SBATCH --error=problems/BipedRobot/4D_model_underactuated/slurm_array/logs/merge_%j.err

# Gathers the chunks into one file. Needs DIONYSOS_NCHUNKS (exported by `submit_all.sh`, which sbatch passes on). The memory is
# about 2 GB + 300 bytes per transition (the transitions of all the chunks are in memory).

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
: "${DIONYSOS_NCHUNKS:?set DIONYSOS_NCHUNKS to the size of the array}"

echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "DIONYSOS_NCHUNKS=$DIONYSOS_NCHUNKS"
echo "DIONYSOS_TRANSITION_OUTDIR=$DIONYSOS_TRANSITION_OUTDIR"
echo "DIONYSOS_ABSTRACTION_FILE=$DIONYSOS_ABSTRACTION_FILE"

srun --export=ALL julia --project=problems/BipedRobot "$DIR/merge_chunks.jl"
