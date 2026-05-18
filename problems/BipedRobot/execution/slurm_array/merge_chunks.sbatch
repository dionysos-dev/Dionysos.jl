#!/bin/bash
#SBATCH --job-name=dionysos-robot-merge
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=../../out/slurm_logs/merge_%j.out
#SBATCH --error=../../out/slurm_logs/merge_%j.err

set -euo pipefail

cd "$(dirname "$0")/../.."

mkdir -p out/slurm_logs

export DIONYSOS_ROBOT_MODEL="${DIONYSOS_ROBOT_MODEL:-6D}"
export DIONYSOS_NCHUNKS="${DIONYSOS_NCHUNKS:-20}"
export DIONYSOS_TRANSITION_OUTDIR="$PWD/out/$DIONYSOS_ROBOT_MODEL/transitions"
export DIONYSOS_ABSTRACTION_FILE="$PWD/out/$DIONYSOS_ROBOT_MODEL/robot_abstraction.jld2"
export DIONYSOS_SIMPLIFY="${DIONYSOS_SIMPLIFY:-3.0}"
export DIONYSOS_TSTEP="${DIONYSOS_TSTEP:-0.1}"

julia --project=. execution/slurm_array/merge_chunks.jl