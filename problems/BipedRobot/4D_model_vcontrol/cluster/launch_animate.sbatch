#!/bin/bash
#SBATCH --job-name=dionysos-robot-animate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=00:20:00
#SBATCH --output=4D_model_vcontrol/cluster/logs/slurm_array/animate_%j.out
#SBATCH --error=4D_model_vcontrol/cluster/logs/slurm_array/animate_%j.err
 
set -euo pipefail
 
mkdir -p 4D_model_vcontrol/cluster/logs
 
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

export DIONYSOS_ABSTRACTION_FILE="$GLOBALSCRATCH/optimizer.jld2"
 
echo "PWD=$(pwd)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "DIONYSOS_ABSTRACTION_FILE=$DIONYSOS_ABSTRACTION_FILE"
 
srun --export=ALL julia --project=. 4D_model_vcontrol/robot_simulation.jl