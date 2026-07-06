#!/bin/bash
#SBATCH --job-name=dionysos-robot-chunks
#SBATCH --array=1-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=4D_model_vcontrol/cluster/logs/slurm_array/chunk_%A_%a.out
#SBATCH --error=4D_model_vcontrol/cluster/logs/slurm_array/chunk_%A_%a.err
 
set -euo pipefail

mkdir -p 4D_model_vcontrol/cluster/logs/
 
mkdir -p 4D_model_vcontrol/cluster/logs/
mkdir -p 4D_model_vcontrol/cluster/logs/transitions
 
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
 
export DIONYSOS_TRANSITION_OUTDIR="4D_model_vcontrol/cluster/logs/transitions"
 
echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT"
echo "DIONYSOS_TRANSITION_OUTDIR=$DIONYSOS_TRANSITION_OUTDIR"
 
srun --export=ALL julia --project=. 4D_model_vcontrol/cluster/compute_chunks.jl