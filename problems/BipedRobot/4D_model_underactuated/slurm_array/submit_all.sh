#!/bin/bash
# Submits the array of chunks and, when all of them have succeeded, the merge:
#
#   problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh NCHUNKS [MAX_PARALLEL]
#
# from the root of the repository. NCHUNKS is the single place where the number of chunks is chosen
# (the array size and the merge both use it); MAX_PARALLEL limits the tasks running at once.
# DIONYSOS_TRANSITION_OUTDIR and DIONYSOS_ABSTRACTION_FILE (environment) choose the directories.

set -euo pipefail

NCHUNKS="${1:?usage: submit_all.sh NCHUNKS [MAX_PARALLEL]}"
THROTTLE="${2:+%$2}"
DIR=problems/BipedRobot/4D_model_underactuated/slurm_array
[ -d "$DIR" ] || { echo "run this from the root of the repository" >&2; exit 1; }

export DIONYSOS_NCHUNKS="$NCHUNKS"
export DIONYSOS_TRANSITION_OUTDIR="${DIONYSOS_TRANSITION_OUTDIR:-$PWD/$DIR/out/transitions}"
mkdir -p "$DIR/logs" "$DIONYSOS_TRANSITION_OUTDIR"

# Chunks left by an earlier run would be merged with the new ones without any warning.
if compgen -G "$DIONYSOS_TRANSITION_OUTDIR/transitions_part_*" > /dev/null; then
    echo "$DIONYSOS_TRANSITION_OUTDIR already has chunk files: remove them (or choose another" >&2
    echo "DIONYSOS_TRANSITION_OUTDIR) so that they are not mixed with this run." >&2
    exit 1
fi

CHUNKS_JOB=$(sbatch --parsable --array="1-${NCHUNKS}${THROTTLE}" "$DIR/submit_compute_chunks.sbatch")
echo "chunks: job $CHUNKS_JOB (array 1-${NCHUNKS}${THROTTLE})"
MERGE_JOB=$(sbatch --parsable --dependency="afterok:${CHUNKS_JOB}" \
    "$DIR/submit_merge_chunks.sbatch")
echo "merge:  job $MERGE_JOB (after $CHUNKS_JOB)"
