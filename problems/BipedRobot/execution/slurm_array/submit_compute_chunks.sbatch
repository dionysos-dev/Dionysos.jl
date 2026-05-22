#!/bin/bash
#SBATCH --job-name=dionysos-robot-chunks
#SBATCH --array=1-200
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm_array/chunk_%A_%a.out
#SBATCH --error=logs/slurm_array/chunk_%A_%a.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs/slurm_array

cd problems/BipedRobot

mkdir -p out/slurm_logs
mkdir -p "out/${DIONYSOS_ROBOT_MODEL:-6D}/transitions"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

export DIONYSOS_ROBOT_MODEL="${DIONYSOS_ROBOT_MODEL:-6D}"
export DIONYSOS_TRANSITION_OUTDIR="$PWD/out/$DIONYSOS_ROBOT_MODEL/transitions"

echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT"
echo "DIONYSOS_TRANSITION_OUTDIR=$DIONYSOS_TRANSITION_OUTDIR"

srun --export=ALL julia --project=. execution/slurm_array/compute_chunk.jl