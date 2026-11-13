#!/usr/bin/env bash
# Runs `ground_cells_linearized.jl` on the configurations of the ground-cells presentation:
# two cell sizes (0.1 and 0.05 rad) times six state boxes (the default box of the script, the cubes
# [-1.5, 1.5]^4 and [-2, 2]^4, and three "knee boxes" with the knees bending backwards,
# q2 in [0, U], q4 in [-U, 0], q1 and q3 in [-U, U], U = 1, 1.5, 2), with five ranges of the foot
# height (the current rule's worst-case margin, the first-order linearization, the linearization plus
# the worst-case second-order remainder, the same with the quadratic term of every cell instead
# (`second_order_cell`), and the corner values).
# Each run writes a plot, a log and a results file (for `compare_ground_cells.jl`) in
# `ground_cells_plots/`, next to this script, and a run whose results file exists is skipped.
#
#   JULIA_ENV=<env> tests/run_ground_cells_matrix.sh [config ...]
#
# `JULIA_ENV` is a Julia environment with Dionysos, Symbolics, Plots, LazySets, StaticArrays and
# MathOptInterface (the test environment by default). With no argument all the configurations run;
# otherwise only the listed ones (`c010_box c010_b15 c010_b20 c005_box c005_b15 c005_b20`, and
# `c010_k10 c010_k15 c010_k20 c005_k10 c005_k15 c005_k20` for the knee boxes). Set
# `METHODS="first_order second_order_cell"` to run only some methods, `FORCE=1` to recompute
# existing results.
# The big boxes take a few minutes per run: their pixels are 1 cm and their sweep steps 0.002–0.0025 rad.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
JULIA_ENV="${JULIA_ENV:-$HERE/../../../../test}"
OUT="$HERE/ground_cells_plots"
mkdir -p "$OUT/logs" "$OUT/wrappers"

# name  cell_size  bounds  pixel_size  sweep_step     (bounds: `box` = the default of the script, a number
# `B` = the cube [-B, B]^4, `kneeU` = the knee box of size U: q2 in [0, U], q4 in [-U, 0], q1, q3 in [-U, U])
CONFIGS=(
    "c010_box 0.1 box 0.005 0.002"
    "c005_box 0.05 box 0.005 0.002"
    "c010_b15 0.1 1.5 0.01 0.002"
    "c005_b15 0.05 1.5 0.01 0.002"
    "c010_b20 0.1 2.0 0.01 0.0025"
    "c005_b20 0.05 2.0 0.01 0.0025"
    "c010_k10 0.1 knee1.0 0.005 0.002"
    "c005_k10 0.05 knee1.0 0.005 0.002"
    "c010_k15 0.1 knee1.5 0.01 0.002"
    "c005_k15 0.05 knee1.5 0.01 0.002"
    "c010_k20 0.1 knee2.0 0.01 0.0025"
    "c005_k20 0.05 knee2.0 0.01 0.0025"
)
# method label : range method : second order
METHOD_SPECS=(
    "margin :margin false"
    "first_order :linearized false"
    "second_order :linearized true"
    "second_order_cell :linearized :cell"
    "corners :corners false"
)

WANTED=("$@")
want_config() {
    [ ${#WANTED[@]} -eq 0 ] && return 0
    for w in "${WANTED[@]}"; do [ "$w" = "$1" ] && return 0; done
    return 1
}
want_method() {
    [ -z "${METHODS:-}" ] && return 0
    for m in $METHODS; do [ "$m" = "$1" ] && return 0; done
    return 1
}

for config in "${CONFIGS[@]}"; do
    read -r name cell bounds pixel step <<<"$config"
    want_config "$name" || continue
    for spec in "${METHOD_SPECS[@]}"; do
        read -r label method second <<<"$spec"
        want_method "$label" || continue
        run="${name}_${label}"
        if [ -f "$OUT/$run.jls" ] && [ -z "${FORCE:-}" ]; then
            echo "== $run: results exist, skipped"
            continue
        fi
        wrapper="$OUT/wrappers/$run.jl"
        {
            echo "using StaticArrays"
            echo "const CELL_SIZE = SVector($cell, $cell, $cell, $cell)"
            if [[ "$bounds" == knee* ]]; then
                u="${bounds#knee}"
                echo "const STATE_BOUNDS = (SVector(-$u, 0.0, -$u, -$u), SVector($u, $u, $u, 0.0))"
            elif [ "$bounds" != "box" ]; then
                echo "const STATE_BOUNDS = (SVector(-$bounds, -$bounds, -$bounds, -$bounds), SVector($bounds, $bounds, $bounds, $bounds))"
            fi
            echo "const RANGE_METHOD = $method"
            echo "const SECOND_ORDER = $second"
            echo "const PIXEL_SIZE = $pixel"
            echo "const SWEEP_STEP = $step"
            echo "const OUTPUT_FILE = \"$OUT/$run.png\""
            echo "const RESULTS_FILE = \"$OUT/$run.jls\""
            echo "include(\"$HERE/ground_cells_linearized.jl\")"
        } >"$wrapper"
        echo "== $run ($(date +%H:%M:%S))"
        julia -t auto --project="$JULIA_ENV" "$wrapper" >"$OUT/logs/$run.log" 2>&1 ||
            echo "   FAILED, see $OUT/logs/$run.log"
        grep -E "^  total" "$OUT/logs/$run.log" || true
    done
done
