# ==============================================================================
# Control synthesis and animation of the normal-form biped experiment, from the saved abstraction
# ==============================================================================
#
# Runs `examples/BipedRobot/biped_4d_normal_form.jl` on the abstraction file written by
# `merge_chunks.jl` (`ABSTRACTION_FILE`): nothing is computed but the synthesis, the closed loop and the
# figures, which are written to `DIONYSOS_OUTPUT_DIR`. The file is never recomputed here: if it is
# missing, this stops (it would be a sequential computation of the whole abstraction).
#
# Rerun it as often as needed with other starts or tolerances: they do not enter the abstraction.
# Put them in a Julia file and point `DIONYSOS_SOLVE_PARAMS` at it; it is included before the example,
# so it sets the variables the example reads, e.g.
#
#     FOOTHOLD_TOLERANCE = 0.03
#     INITIAL_VELOCITY = (-0.44, -0.36)       # default: any start momentum
#     Q0 = SVector(0.0, 0.5, -0.4, -0.9)      # StaticArrays is loaded for it
#     TARGET_VELOCITY = (-0.3, -0.1)          # default: any end momentum
#     START_FOOTHOLDS = [SVector(x, 0.0) for x in -0.25:0.01:-0.15]   # start: foot on these points
#
# (`FOOTHOLD_TOLERANCE`, `Q0`, `INITIAL_VELOCITY`, `TARGET_VELOCITY`, `START_FOOTHOLDS`.) What defines
# the abstraction (`DISC`, `X_BOUNDS`, `ONE_JOINT_PER_STEP`, `OBSTACLE`) is the example's: the file is
# refused if it differs.
#
# The start is the configuration `Q0` with any momentum of `INITIAL_VELOCITY`, or, with
# `START_FOOTHOLDS`, every configuration with the swing foot at one of these points: the solve succeeds
# if one of the start cells is in the controller's domain. The closed loop runs from each of them, and
# one is animated (the fewest worst-case steps; with `START_FOOTHOLDS`, the foot closest to their
# middle first). The solver's own `success` asks for all of them.
#
# Environment:
#   DIONYSOS_ABSTRACTION_FILE   the saved abstraction (default: slurm_array/out/normal_form_abstraction.jld2)
#   DIONYSOS_PARAMS             the parameter file of the chunks, if any (see compute_chunk.jl), included
#                               first
#   DIONYSOS_SOLVE_PARAMS       optional Julia file with the variables above
#   DIONYSOS_SCENARIO           name of the outputs (default: job<SLURM_JOB_ID>)
#   DIONYSOS_OUTPUT_DIR         directory of the outputs (default: slurm_array/out/results)

# `instantiate` must not precompile the whole project (MeshCat, RigidBodyDynamics, ... that the experiment
# never loads): the precompilation is the job of `precompile.jl`, once, before the array.
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
import Pkg
Pkg.instantiate()

using StaticArrays   # for the params file
using Dates
import JLD2

stamp(msg) = @info "[solve $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] $msg"
peak_gb() = round(Sys.maxrss() / 2^30; digits = 1)

# The data of the figures of a report, saved by every solve, with or without a controller (the globals
# are those of the example): the numbers of `RP.synthesis_report`, the cells of the domain with their
# worst-case number of steps, the start cells (the momentum cells of the start configuration, with
# their worst-case steps, `-1` outside the domain) and the closed loops from those in the domain, the
# optimistic path from `x_start` (the start of the animated run; cell centres and inputs) and its
# replay on the RK4 map with the same inputs, and the closed loops from the deepest cell of the domain
# (most steps to the target, the swing foot farthest back among them) and from the cell with the swing
# foot farthest back. Plotted by `tests/plot_step_diagnostics.jl`.
# The synthesized controller, to rerun closed loops without solving again: for every cell of the
# domain outside the target, its grid position (`state_origin + state_h .* position` is the centre)
# and its input as multiples of `du` (`input_steps .* du`). About 13 bytes per cell.
function controller_table()
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    positions = NTuple{5, Int16}[]
    input_steps = NTuple{3, Int8}[]
    for q in controller.dom
        us = controller.controller_map(q)
        isempty(us) && continue          # a target cell: no input
        x = SY.get_concrete_state(abstract_system, q)
        u = SY.get_concrete_input(abstract_system, first(us))
        push!(positions, Int16.(Tuple(MP.get_pos_by_coord(disc.state_grid, x))))
        push!(input_steps, Int8.(Tuple(round.(Int, u ./ disc.du))))
    end
    return (;
        positions,
        input_steps,
        state_origin = SVector{5}(MP.get_origin(disc.state_grid)),
        state_h = SVector{5}(MP.get_h(disc.state_grid)),
        disc.du,
    )
end

function save_diagnostics(file)
    as_run(start, steps, traj) = (;
        start,
        steps,
        states = collect(ST.states(traj)),
        inputs = ST.inputs(traj) === nothing ? SVector{3, Float64}[] :
                 collect(ST.inputs(traj)),
    )
    path = report.optimistic_path
    replay = nothing
    if path !== nothing
        replay = [SVector{5}(x_start)]
        for u in path.inputs
            push!(replay, concrete_system.f(replay[end], u))
        end
    end
    foot_x(c) = RP.swing_foot_position(SVector{4}(c[2:5]), geometry)[1]
    outside = [c for c in report.domain_cells if c[2] > 0]
    runs = Dict{String, Any}()
    if !isempty(outside)
        deepest_steps = maximum(last, outside)
        # (Not `starts`: that is the example's start cells, read below.)
        extremes = (
            "deepest" => argmin(
                c -> foot_x(c[1]),
                [c for c in outside if c[2] == deepest_steps],
            ),
            "farthest_back" => argmin(c -> foot_x(c[1]), outside),
        )
        for (name, (start, steps)) in extremes
            # The same cell twice (`Base.`: `values` is the example's parameters).
            any(r -> r.start == start, Base.values(runs)) && continue
            traj = RP.closed_loop(optimizer, step_pb, start; ankle_name = "p₁")
            runs[name] = as_run(start, steps, traj)
        end
    end
    JLD2.jldsave(
        file;
        counts = (;
            report.target,
            report.carved_targets,
            report.domain,
            report.optimistic,
            report.x0_in_domain,
            report.x0_in_optimistic,
        ),
        domain_steps = report.domain_steps,
        optimistic_steps = report.optimistic_steps,
        domain_cells = report.domain_cells,
        optimistic_path = path,
        replay,
        runs,
        start_cells = [(SVector{5}(c.state), c.steps) for c in starts.cells],
        start_runs = [as_run(SVector{5}(r.state), r.steps, r.traj) for r in starts.runs],
        controller = controller_table(),
        x0 = SVector{5}(x_start),   # the start of the report and of the animated run
        x0_reference = SVector{5}(x0),
        tstep = disc.tstep,
        x_lb = SVector{5}(x_lb),
        x_ub = SVector{5}(x_ub),
        foothold_tolerance,
    )
    return file
end

const EXAMPLE = abspath(
    joinpath(
        @__DIR__,
        "..",
        "..",
        "..",
        "..",
        "examples",
        "BipedRobot",
        "biped_4d_normal_form.jl",
    ),
)

try
    stamp(
        "start: job $(get(ENV, "SLURM_JOB_ID", "-")), memory requested " *
        "$(get(ENV, "SLURM_MEM_PER_NODE", "?")) MB, Julia $(VERSION), $(gethostname())",
    )
    abstraction_path = abspath(
        get(
            ENV,
            "DIONYSOS_ABSTRACTION_FILE",
            joinpath(@__DIR__, "out", "normal_form_abstraction.jld2"),
        ),
    )
    isfile(abstraction_path) || error(
        "$abstraction_path does not exist: run the chunks and the merge first (submit_all.sh), " *
        "or set DIONYSOS_ABSTRACTION_FILE",
    )
    n = JLD2.jldopen(
        f -> haskey(f, "n_transitions") ? f["n_transitions"] : nothing,
        abstraction_path,
        "r",
    )
    stamp(
        "abstraction $abstraction_path: $(round(filesize(abstraction_path) / 2^30; digits = 1)) GB, " *
        "$(n === nothing ? "?" : n) transitions (in memory, an automaton takes of the order of " *
        "$(n === nothing ? "?" : round(Int, 200 * n / 2^30)) GB for them)",
    )

    # Read by the example before it defaults them.
    global ABSTRACTION_FILE = abstraction_path
    global SCENARIO_NAME =
        Symbol(get(ENV, "DIONYSOS_SCENARIO", "job" * get(ENV, "SLURM_JOB_ID", "local")))
    global OUTPUT_DIR =
        abspath(get(ENV, "DIONYSOS_OUTPUT_DIR", joinpath(@__DIR__, "out", "results")))
    global RUN_PLOTS = true

    # The parameters of the experiment (those of the chunks), then the ones of this solve, which may
    # override them.
    experiment = get(ENV, "DIONYSOS_PARAMS", "")
    if !isempty(experiment)
        isfile(experiment) || error("DIONYSOS_PARAMS = $experiment does not exist")
        include(abspath(experiment))
    end
    params = get(ENV, "DIONYSOS_SOLVE_PARAMS", "")
    if !isempty(params)
        isfile(params) || error("DIONYSOS_SOLVE_PARAMS = $params does not exist")
        include(abspath(params))      # may override the variables above too (e.g. SCENARIO_NAME)
    end

    mkpath(OUTPUT_DIR)
    prefix = joinpath(OUTPUT_DIR, "biped_4d_normal_form_$(SCENARIO_NAME)")
    # A record of what produced it.
    isempty(experiment) || cp(experiment, prefix * "_experiment.jl"; force = true)
    isempty(params) || cp(params, prefix * "_params.jl"; force = true)

    stamp(
        "setting up the model, then loading the abstraction; scenario $SCENARIO_NAME, " *
        "outputs $prefix*, params '$params'",
    )
    t = @elapsed include(EXAMPLE)
    stamp(
        "synthesis and closed loop done in $(round(t; digits = 1)) s (peak memory $(peak_gb()) GB)",
    )
    try
        stamp("diagnostics saved to $(save_diagnostics(prefix * "_diagnostics.jld2"))")
    catch e
        @warn "the diagnostics were not saved" exception = (e, catch_backtrace())
    end

    if isfile(prefix * ".gif")
        reached = count(r -> last(ST.states(r.traj)) ∈ step_pb.target_set, starts.runs)
        stamp(
            "SOLVE OK: $(length(starts.runs)) of $(length(starts.cells)) start cells in the " *
            "domain, $reached closed loops reach the target; the run from " *
            "p₁ = $(round(x_start[1]; digits = 4)): $(prefix * ".gif")",
        )
    else
        stamp(
            "SOLVE FAILED: no animation was written: " * (
                isempty(starts.runs) ?
                "no start cell (Q0 or START_FOOTHOLDS, with the momenta of INITIAL_VELOCITY) is " *
                "in the controller's domain. Try a larger FOOTHOLD_TOLERANCE or another start." :
                "the run from the start has no step (its cell is in the target already, or " *
                "the closed loop did not run: see the lines above)."
            ),
        )
        exit(2)
    end
catch e
    stamp(
        "SOLVE FAILED after peak memory $(peak_gb()) GB: " * (
            e isa OutOfMemoryError ?
            "OutOfMemoryError: raise --mem of submit_solve.sbatch" : sprint(showerror, e)
        ),
    )
    @error "SOLVE FAILED" exception = (e, catch_backtrace())
    exit(1)
end
