# Sourced by the sbatch scripts of the pipeline (`source "$DIR/env.sh"`): where the Julia precompilation
# cache of the pipeline lives.
#
# The cache holds machine code, built for the CPU that compiled it, and a cache file has the same name whatever
# the CPU: on a cluster whose nodes are not all alike and share the home directory, jobs on different CPU types
# overwrite each other's cache files, each one finding the cache of the other invalid, and all of them
# recompile, waiting for each other's lock files ("Being precompiled by another machine"). So the pipeline
# keeps its own depot, in front of ~/.julia (which is still where the sources of the packages are found), with
# one cache directory per kind of CPU: jobs on the same CPU type share it, the others have their own.
#   - With JULIA_CPU_TARGET set (a portable target, see README §1), all the nodes share one directory.
#   - DIONYSOS_CACHE_DEPOT: where these directories are (default: ~/.julia_biped_cache).

# sbatch passes on the environment of the shell that submitted the job, PATH included: a `julia` that only an
# interactive shell can find (juliaup adds it in ~/.bash_profile) is not found by the job.
command -v julia > /dev/null 2>&1 || {
    echo "julia is not in the PATH of this job. Submit from a shell where 'julia' works (a login shell: ssh host 'bash -l -c ...'), or put its directory in the PATH before submitting." >&2
    exit 1
}

DIONYSOS_BASE_DEPOT="${DIONYSOS_BASE_DEPOT:-${JULIA_DEPOT_PATH:-$HOME/.julia}}"
if [ -n "${JULIA_CPU_TARGET:-}" ]; then
    DIONYSOS_CACHE_KEY=portable
else
    DIONYSOS_CACHE_KEY="$(julia --startup-file=no -e 'print(Sys.CPU_NAME)')"
fi
export DIONYSOS_BASE_DEPOT DIONYSOS_CACHE_KEY
export JULIA_DEPOT_PATH="${DIONYSOS_CACHE_DEPOT:-$HOME/.julia_biped_cache}/$DIONYSOS_CACHE_KEY:$DIONYSOS_BASE_DEPOT"
echo "node $(hostname), CPU type $DIONYSOS_CACHE_KEY, JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH"
