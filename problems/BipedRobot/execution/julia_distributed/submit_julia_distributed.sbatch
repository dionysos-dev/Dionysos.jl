#!/bin/bash
#SBATCH --job-name=dionysos_biped_dist
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=64G
#SBATCH --time=00:10:00

set -euo pipefail

mkdir -p logs

export DIONYSOS_ROBOT_MODEL=6D
export DIONYSOS_NWORKERS=11
export DIONYSOS_NPARTS=10
export DIONYSOS_PARTITION_STRATEGY=contiguous
export DIONYSOS_SIMPLIFY=1.0
export DIONYSOS_TSTEP=0.1

echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "DIONYSOS_NWORKERS=$DIONYSOS_NWORKERS"
echo "DIONYSOS_NPARTS=$DIONYSOS_NPARTS"

srun --export=ALL julia --project=problems/BipedRobot \
    problems/BipedRobot/execution/julia_distributed/run_abstraction.jl