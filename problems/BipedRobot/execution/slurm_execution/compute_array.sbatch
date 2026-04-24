#!/bin/bash
#SBATCH --job-name=dionysos-robot-chunks
#SBATCH --array=1-20
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=../../out/slurm_logs/chunk_%A_%a.out
#SBATCH --error=../../out/slurm_logs/chunk_%A_%a.err

set -euo pipefail

cd "$(dirname "$0")/../.."

mkdir -p out/slurm_logs
mkdir -p out/${DIONYSOS_ROBOT_MODEL:-6D}/transitions

export DIONYSOS_ROBOT_MODEL="${DIONYSOS_ROBOT_MODEL:-6D}"
export DIONYSOS_TRANSITION_OUTDIR="$PWD/out/$DIONYSOS_ROBOT_MODEL/transitions"
export DIONYSOS_SIMPLIFY="${DIONYSOS_SIMPLIFY:-3.0}"
export DIONYSOS_TSTEP="${DIONYSOS_TSTEP:-0.1}"

julia --project=. execution/slurm_array/compute_chunk.jl