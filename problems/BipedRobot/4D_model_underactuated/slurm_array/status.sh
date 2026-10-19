#!/bin/bash
# What happened to the jobs of a pipeline: state, exit code, peak memory, time, and the end of the logs.
#
#   problems/BipedRobot/4D_model_underactuated/slurm_array/status.sh            # the last submit_all.sh
#   problems/BipedRobot/4D_model_underactuated/slurm_array/status.sh CHUNKS_JOB MERGE_JOB [SOLVE_JOB]
#
# from the root of the repository. In the State column: COMPLETED is fine; OUT_OF_MEMORY means the job
# needed more than its `--mem` (a job killed by Slurm cannot write anything: its log just stops);
# TIMEOUT more than its `--time`; CANCELLED by a failed dependency is the consequence of an earlier
# failure, look at the first job that is not COMPLETED.

DIR=problems/BipedRobot/4D_model_underactuated/slurm_array
if [ $# -ge 2 ]; then
    CHUNKS_JOB=$1; MERGE_JOB=$2; SOLVE_JOB=${3:-none}
else
    [ -f "$DIR/logs/last_submission.txt" ] || { echo "no $DIR/logs/last_submission.txt: give the job ids" >&2; exit 1; }
    eval "$(awk '{for (i = 1; i <= NF; i++) { split($i, a, "="); \
        if (a[1] == "chunks") print "CHUNKS_JOB=" a[2]; \
        if (a[1] == "merge") print "MERGE_JOB=" a[2]; \
        if (a[1] == "solve") print "SOLVE_JOB=" a[2] }}' "$DIR/logs/last_submission.txt")"
fi
FMT=JobID%24,JobName%20,State%20,ExitCode,Elapsed,ReqMem,MaxRSS

echo "== chunks (job $CHUNKS_JOB): tasks by state"
sacct -X -n -j "$CHUNKS_JOB" -o State%24 | sed 's/ *$//' | sort | uniq -c
FAILED=$(sacct -X -n -j "$CHUNKS_JOB" -o JobID%20,State%20 | awk '$2!="COMPLETED"{print $1}' | head -20 | tr '\n' ' ')
[ -n "$FAILED" ] && echo "tasks not COMPLETED (first 20): $FAILED"
echo "chunk files: $(ls "${DIONYSOS_TRANSITION_OUTDIR:-$DIR/out/transitions}" 2>/dev/null | grep -c '^transitions_part_')"

for JOB in "$MERGE_JOB" "$SOLVE_JOB"; do
    [ "$JOB" = none ] && continue
    echo; echo "== job $JOB"
    sacct -j "$JOB" --format=$FMT
    for LOG in "$DIR"/logs/*_"$JOB".out "$DIR"/logs/*_"$JOB".err; do
        [ -f "$LOG" ] || continue
        echo "-- last lines of $LOG"; tail -n 8 "$LOG"
    done
done
echo
echo "results: $(ls "${DIONYSOS_OUTPUT_DIR:-$DIR/out/results}"/*.gif 2>/dev/null | tr '\n' ' ')"
