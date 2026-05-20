#!/bin/bash
#SBATCH --job-name=dionysos-robot-merge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=problems/BipedRobot/out/slurm_logs/merge_%j.out
#SBATCH --error=problems/BipedRobot/out/slurm_logs/merge_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

mkdir -p problems/BipedRobot/out/slurm_logs

cd problems/BipedRobot

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

export DIONYSOS_ROBOT_MODEL="${DIONYSOS_ROBOT_MODEL:-6D}"
export DIONYSOS_NCHUNKS="${DIONYSOS_NCHUNKS:-200}"
export DIONYSOS_TRANSITION_OUTDIR="$PWD/out/$DIONYSOS_ROBOT_MODEL/transitions"
export DIONYSOS_ABSTRACTION_FILE="$PWD/out/$DIONYSOS_ROBOT_MODEL/robot_abstraction.jld2"
export DIONYSOS_SIMPLIFY="${DIONYSOS_SIMPLIFY:-1.0}"
export DIONYSOS_TSTEP="${DIONYSOS_TSTEP:-0.1}"

echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "DIONYSOS_ROBOT_MODEL=$DIONYSOS_ROBOT_MODEL"
echo "DIONYSOS_NCHUNKS=$DIONYSOS_NCHUNKS"
echo "DIONYSOS_TRANSITION_OUTDIR=$DIONYSOS_TRANSITION_OUTDIR"
echo "DIONYSOS_ABSTRACTION_FILE=$DIONYSOS_ABSTRACTION_FILE"

srun --export=ALL julia --project=. execution/slurm_array/merge_chunks.jl