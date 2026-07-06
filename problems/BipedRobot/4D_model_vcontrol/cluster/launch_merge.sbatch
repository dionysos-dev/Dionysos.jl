#!/bin/bash
#SBATCH --job-name=dionysos-robot-merge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --output=4D_model_vcontrol/cluster/logs/slurm_array/merge_%j.out
#SBATCH --error=4D_model_vcontrol/cluster/logs/slurm_array/merge_%j.err
 
set -euo pipefail
 
mkdir -p 4D_model_vcontrol/cluster/out
 
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
 
export DIONYSOS_TRANSITION_OUTDIR="4D_model_vcontrol/cluster/logs/transitions"
export DIONYSOS_ABSTRACTION_FILE="4D_model_vcontrol/cluster/out/optimizer.jld2"
 
echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "DIONYSOS_TRANSITION_OUTDIR=$DIONYSOS_TRANSITION_OUTDIR"
echo "DIONYSOS_ABSTRACTION_FILE=$DIONYSOS_ABSTRACTION_FILE"
 
srun --export=ALL julia --project=. 4D_model_vcontrol/cluster/merge_chunks.jl