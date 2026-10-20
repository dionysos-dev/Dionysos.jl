#!/bin/bash
# Submits the whole pipeline: the precompilation of the packages on a compute node, then (each after the
# previous one has succeeded) the array of chunks, the merge into one abstraction file, and the control
# synthesis, closed loop and animation:
#
#   problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh NCHUNKS [MAX_PARALLEL]
#
# from the root of the repository. NCHUNKS is the single place where the number of chunks is chosen
# (the array size and the merge both use it); MAX_PARALLEL limits the tasks running at once.
# DIONYSOS_TRANSITION_OUTDIR, DIONYSOS_ABSTRACTION_FILE and DIONYSOS_OUTPUT_DIR (environment) choose the
# directories; DIONYSOS_SOLVE_PARAMS and DIONYSOS_SCENARIO are passed to the last job (see solve.jl).
# With DIONYSOS_SKIP_SOLVE=1 only the chunks and the merge are submitted: the abstraction file is saved,
# and `sbatch submit_solve.sbatch` runs the synthesis (again) later, with other starts or tolerances.

set -euo pipefail

NCHUNKS="${1:?usage: submit_all.sh NCHUNKS [MAX_PARALLEL]}"
command -v julia > /dev/null 2>&1 || { echo "julia is not in the PATH: the jobs inherit it, submit from a shell where 'julia' works (a login shell)" >&2; exit 1; }
THROTTLE="${2:+%$2}"
DIR=problems/BipedRobot/4D_model_underactuated/slurm_array
[ -d "$DIR" ] || { echo "run this from the root of the repository" >&2; exit 1; }

# Extra arguments for every sbatch call, e.g. to keep the whole pipeline on one kind of node (the Julia
# precompilation cache is per CPU type, README §6):
#   DIONYSOS_SBATCH_ARGS="--partition=cp3 --constraint=Skylake" submit_all.sh 300
# (the features of the nodes: sinfo -p <partition> -o "%20N %f"; or --exclude=<nodes>, --nodelist=<nodes>).
read -r -a EXTRA <<< "${DIONYSOS_SBATCH_ARGS:-}"
SBATCH=(sbatch --parsable ${EXTRA[@]+"${EXTRA[@]}"})
[ -n "${DIONYSOS_SBATCH_ARGS:-}" ] && echo "sbatch arguments for every job: $DIONYSOS_SBATCH_ARGS"

# The Julia cache is per CPU type and the precompile job builds it for the node it lands on: if the nodes of the
# partition are not all alike and nothing pins the jobs to one type, the tasks of the other types find no
# cache (and stop). So this is refused unless the jobs are pinned (--constraint, --nodelist, --exclude), or a
# portable JULIA_CPU_TARGET is set, or DIONYSOS_ANY_NODE=1.
if [[ " ${DIONYSOS_SBATCH_ARGS:-} " != *" --constraint"* && " ${DIONYSOS_SBATCH_ARGS:-} " != *" -C"* \
   && " ${DIONYSOS_SBATCH_ARGS:-} " != *" --nodelist"* && " ${DIONYSOS_SBATCH_ARGS:-} " != *" --exclude"* \
   && -z "${JULIA_CPU_TARGET:-}" && "${DIONYSOS_ANY_NODE:-0}" != "1" ]]; then
    NFEATURES=$(sinfo -h -o "%f" 2>/dev/null | sort -u | wc -l)
    if [ "$NFEATURES" -gt 1 ]; then
        echo "The cluster has nodes of different kinds ($NFEATURES feature sets):" >&2
        sinfo -h -o "  %20N %c cores  %f" | sort -u >&2
        echo "and the Julia precompilation cache is per CPU type: pin every job to one of them, e.g." >&2
        echo "  DIONYSOS_SBATCH_ARGS=\"--constraint=<one feature, e.g. Rome>\" $0 $NCHUNKS" >&2
        echo "(or export a portable JULIA_CPU_TARGET, or DIONYSOS_ANY_NODE=1 to submit anyway: README section 1)." >&2
        exit 1
    fi
fi

export DIONYSOS_NCHUNKS="$NCHUNKS"
export DIONYSOS_TRANSITION_OUTDIR="${DIONYSOS_TRANSITION_OUTDIR:-$PWD/$DIR/out/transitions}"
mkdir -p "$DIR/logs" "$DIONYSOS_TRANSITION_OUTDIR"

# Chunks left by an earlier run would be merged with the new ones without any warning.
if compgen -G "$DIONYSOS_TRANSITION_OUTDIR/transitions_part_*" > /dev/null; then
    echo "$DIONYSOS_TRANSITION_OUTDIR already has chunk files: remove them (or choose another" >&2
    echo "DIONYSOS_TRANSITION_OUTDIR) so that they are not mixed with this run." >&2
    exit 1
fi

# The precompilation cache is built once, on a compute node, before the array: otherwise every task
# compiles the packages itself (DIONYSOS_SKIP_PRECOMPILE=1 to skip it, e.g. if it is known to be valid).
PRE_JOB=none
PRE_DEP=()
if [ "${DIONYSOS_SKIP_PRECOMPILE:-0}" != "1" ]; then
    PRE_JOB=$("${SBATCH[@]}" "$DIR/submit_precompile.sbatch")
    echo "precompile: job $PRE_JOB"
    PRE_DEP=(--dependency="afterok:${PRE_JOB}" --kill-on-invalid-dep=yes)
fi
CHUNKS_JOB=$("${SBATCH[@]}" ${PRE_DEP[@]+"${PRE_DEP[@]}"} --array="1-${NCHUNKS}${THROTTLE}" "$DIR/submit_compute_chunks.sbatch")
echo "chunks: job $CHUNKS_JOB (array 1-${NCHUNKS}${THROTTLE})"
# kill-on-invalid-dep: if a job fails, the ones that wait for it are cancelled instead of staying pending.
MERGE_JOB=$("${SBATCH[@]}" --dependency="afterok:${CHUNKS_JOB}" --kill-on-invalid-dep=yes \
    "$DIR/submit_merge_chunks.sbatch")
echo "merge:  job $MERGE_JOB (after $CHUNKS_JOB)"
SOLVE_JOB=none
if [ "${DIONYSOS_SKIP_SOLVE:-0}" = "1" ]; then
    echo "solve:  skipped (DIONYSOS_SKIP_SOLVE=1)"
else
    SOLVE_JOB=$("${SBATCH[@]}" --dependency="afterok:${MERGE_JOB}" --kill-on-invalid-dep=yes \
        "$DIR/submit_solve.sbatch")
    echo "solve:  job $SOLVE_JOB (after $MERGE_JOB)"
fi
echo "precompile=$PRE_JOB chunks=$CHUNKS_JOB merge=$MERGE_JOB solve=$SOLVE_JOB nchunks=$NCHUNKS" > "$DIR/logs/last_submission.txt"
echo "check what happened with: $DIR/status.sh"
[ "$SOLVE_JOB" = none ] && exit 0
echo "results: ${DIONYSOS_OUTPUT_DIR:-$PWD/$DIR/out/results}/biped_4d_normal_form_job${SOLVE_JOB}.gif"
