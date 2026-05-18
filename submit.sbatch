#!/bin/bash
#SBATCH --job-name=dionysos_biped
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=09:00:00

set -euo pipefail

mkdir -p logs
source problems/BipedRobot/setup_julia.sh

echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "DIONYSOS_NWORKERS=$DIONYSOS_NWORKERS"

srun --export=ALL julia --project=problems/BipedRobot problems/BipedRobot/6D_model/robot_example.jl