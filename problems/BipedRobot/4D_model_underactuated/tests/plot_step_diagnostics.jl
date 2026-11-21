# Figures of a footstep run from the diagnostics file that `slurm_array/solve.jl` saves with every solve
# (`<prefix>_diagnostics.jld2`), whether or not a controller was found:
#
# * `<name>_domain.png`: the controller's domain, (a) as the swing foot at the cell centres, coloured by
#   the worst-case number of steps to the target, with the ground, the obstacle, the foothold and the
#   start, and (b) in the ankle plane `(q₁, p₁)`: the cells of the target, those of the rest of the
#   domain, the start cells (the momentum cells of the start configuration, in the domain or not), the
#   start, the optimistic path from the start and its replay on the RK4 map;
# * `<name>_steps.png`: the number of cells per number of steps to the target, worst case (the domain)
#   and optimistic (some successors);
# * `<name>_starts.png`: the closed loops from the start cells in the domain, side by side (swing
#   foot, ankle angle and momentum);
# * `<name>_<run>.png`, `_<run>_ankle.png`, `_<run>.gif` (`RP.plot_footstep`) for the closed loops of
#   the diagnostics file (`deepest`: from the domain cell with the most steps to the target;
#   `farthest_back`: from the domain cell with the swing foot farthest back);
# * `<name>_optimistic.png` and `.gif`: the optimistic path from the start, cell by cell (left), next
#   to the RK4 map from the start with the same inputs (right), and the ankle along both.
#
# The experiment file (the `DIONYSOS_PARAMS` of the run) gives the grid, the box and the obstacle:
#
#   julia --project=problems/BipedRobot problems/BipedRobot/4D_model_underactuated/tests/plot_step_diagnostics.jl \
#       DIAGNOSTICS.jld2 EXPERIMENT.jl [OUTDIR] [NAME]

using Printf
using StaticArrays
import Dionysos
import JLD2
import LazySets
import Plots

const UT = Dionysos.Utils
const MP = Dionysos.Mapping
const ST = Dionysos.System
const PROBLEM_DIR = dirname(@__DIR__)
include(joinpath(PROBLEM_DIR, "robot_problem.jl"))
import .UnderactuatedRobotProblem as RP
include(joinpath(PROBLEM_DIR, "BipedAbstraction.jl"))
import .BipedAbstraction as BA

length(ARGS) >= 2 ||
    error("usage: plot_step_diagnostics.jl DIAGNOSTICS.jld2 EXPERIMENT.jl [OUTDIR] [NAME]")
const DIAG_FILE = ARGS[1]
const OUTDIR = length(ARGS) >= 3 ? ARGS[3] : dirname(abspath(DIAG_FILE))
const NAME =
    length(ARGS) >= 4 ? ARGS[4] : replace(basename(DIAG_FILE), "_diagnostics.jld2" => "")
include(abspath(ARGS[2]))      # DISC_KWARGS, X_BOUNDS, OBSTACLE, GROUND_AMBIGUOUS
mkpath(OUTDIR)
out(suffix) = joinpath(OUTDIR, NAME * suffix)

d = JLD2.load(DIAG_FILE)
x0 = d["x0"]
foothold = SVector(0.2, 0.0)
obstacle = @isdefined(OBSTACLE) ? OBSTACLE : nothing

# The model, the grid, the carving and the dynamics, as in the example.
values_ = RP.default_parameters()
model = RP.symbolic_model()
geometry = RP.Geometry(model, values_)
reduced = RP.reduced_dynamics(model)
f_sym = RP.substitute_parameters(model, reduced.f, values_)
disc = RP.default_discretization(; DISC_KWARGS...)
x_lb, x_ub = X_BOUNDS
X_box = LazySets.Hyperrectangle(; low = x_lb, high = x_ub)
carving = RP.infeasible_cells(
    X_box,
    disc.state_grid,
    obstacle,
    geometry;
    obstacle_name = "lift constraint",
    ground_ambiguous = @isdefined(GROUND_AMBIGUOUS) ? GROUND_AMBIGUOUS : :keep,
)
U_box = LazySets.Hyperrectangle(;
    low = SVector(-disc.u_max, -disc.u_max, -disc.u_max),
    high = SVector(disc.u_max, disc.u_max, disc.u_max),
)
certified = BA.CertifiedRK4GrowthBound(
    f_sym,
    reduced.x,
    reduced.u,
    UT.set_minus(X_box, carving.cells),
    U_box,
    disc.tstep,
)
dynamics = (x, u) -> SVector{5}(certified.f(Vector(x), Vector(u)))
foot(x) = RP.swing_foot_position(SVector{4}(x[2:5]), geometry)
h = MP.get_h(disc.state_grid)
tlabel = @sprintf("tstep = 1/%d s", round(Int, 1 / disc.tstep))

function draw_scene!(fig; xlims = (-0.6, 0.6), ylims = (-0.2, 0.55))
    Plots.plot!(fig, [xlims[1], xlims[2]], [0.0, 0.0]; color = :black, lw = 1, label = "")
    for ob in RP._obstacle_list(obstacle)
        Plots.plot!(fig, ob; color = :gray, alpha = 0.5, label = "")
    end
    Plots.scatter!(
        fig,
        [foothold[1]],
        [foothold[2]];
        marker = :xcross,
        ms = 7,
        color = :green,
        label = "",
    )
    Plots.plot!(fig; xlims, ylims, aspect_ratio = :equal)
    return fig
end

# ------------------------------------------------------------
# The domain: swing foot and ankle plane
# ------------------------------------------------------------
cells = d["domain_cells"]
path = d["optimistic_path"]
replay = d["replay"]
println(
    "$(length(cells)) domain cells; x0 in the domain: $(d["counts"].x0_in_domain), in the optimistic set: " *
    "$(d["counts"].x0_in_optimistic)",
)

# (a) The swing foot depends on the angles only: one point per angle cell, with the most steps of its p₁ copies.
by_angles = Dict{NTuple{4, Int}, Tuple{SVector{2, Float64}, Int}}()
for (c, k) in cells
    key = Tuple(round.(Int, c[2:5] ./ h[2:5]))
    old = get(by_angles, key, (foot(c), -1))
    by_angles[key] = (old[1], max(old[2], k))
end
pts = collect(values(by_angles))
order = sortperm(pts; by = last)
pa = Plots.plot(;
    title = "(a) swing foot at the domain cells",
    xlabel = "x (m)",
    ylabel = "z (m)",
    legend = :topleft,
)
draw_scene!(pa; xlims = (-0.35, 0.5), ylims = (-0.1, 0.3))
Plots.scatter!(
    pa,
    [pts[i][1][1] for i in order],
    [pts[i][1][2] for i in order];
    zcolor = [pts[i][2] for i in order],
    color = :viridis,
    ms = 3,
    msw = 0,
    label = "domain (colour: worst-case steps)",
)
f0 = foot(x0)
Plots.scatter!(
    pa,
    [f0[1]],
    [f0[2]];
    marker = :star5,
    ms = 9,
    color = :red,
    label = "start x0",
)

# (b) The ankle plane: cells of the target and of the rest of the domain, the start, the paths.
pb = Plots.plot(;
    title = "(b) ankle plane",
    xlabel = "q₁ (rad)",
    ylabel = "p₁ (N m s)",
    legend = :right,
)
ankle(c) = (round(c[2] / h[2]) * h[2], round(c[1] / h[1]) * h[1])
target_ankle = unique(ankle(c) for (c, k) in cells if k == 0)
rest_ankle = setdiff(unique(ankle(c) for (c, k) in cells if k > 0), target_ankle)
Plots.scatter!(
    pb,
    first.(target_ankle),
    last.(target_ankle);
    marker = :square,
    ms = 7,
    msw = 0,
    color = :lightgreen,
    label = "target",
)
if isempty(rest_ankle)
    # Every cell of the domain has the (q₁, p₁) cell of a target cell: the ankle never left the target's.
    Plots.annotate!(
        pb,
        0.7 * x_ub[2],
        x_lb[1] + 0.06,
        Plots.text("rest of the domain: no ankle\ncell outside the target's", 8, :blue),
    )
else
    Plots.scatter!(
        pb,
        first.(rest_ankle),
        last.(rest_ankle);
        ms = 3,
        color = :blue,
        label = "rest of the domain",
    )
end
if path !== nothing
    Plots.plot!(
        pb,
        [s[2] for s in path.states],
        [s[1] for s in path.states];
        marker = :circle,
        color = :orange,
        lw = 2,
        label = "optimistic path (cells)",
    )
    Plots.plot!(
        pb,
        [s[2] for s in replay],
        [s[1] for s in replay];
        marker = :diamond,
        color = :purple,
        ls = :dash,
        lw = 2,
        label = "RK4, same inputs",
    )
end
start_cells = get(d, "start_cells", nothing)   # saved since the start is a configuration
if start_cells !== nothing
    out_cells = [c for (c, k) in start_cells if k < 0]
    in_cells = [c for (c, k) in start_cells if k >= 0]
    isempty(out_cells) || Plots.scatter!(
        pb,
        [c[2] for c in out_cells],
        [c[1] for c in out_cells];
        ms = 3,
        color = :white,
        msc = :gray,
        label = "start cells outside the domain",
    )
    isempty(in_cells) || Plots.scatter!(
        pb,
        [c[2] for c in in_cells],
        [c[1] for c in in_cells];
        ms = 4,
        color = :red,
        label = "start cells in the domain",
    )
end
Plots.scatter!(
    pb,
    [x0[2]],
    [x0[1]];
    marker = :star5,
    ms = 10,
    color = :red,
    label = "start x0",
)
Plots.plot!(pb; xlims = (x_lb[2], x_ub[2]), ylims = (x_lb[1], x_ub[1]))
Plots.savefig(
    Plots.plot(
        pa,
        pb;
        layout = Plots.grid(1, 2; widths = [0.56, 0.44]),
        size = (1500, 560),
        left_margin = 6Plots.mm,
        bottom_margin = 5Plots.mm,
        plot_title = "Controller's domain ($tlabel)",
    ),
    out("_domain.png"),
)

# ------------------------------------------------------------
# Number of cells per number of steps
# ------------------------------------------------------------
ds, os = d["domain_steps"], d["optimistic_steps"]
K = max(length(ds), length(os))
pad(v) = vcat(v, zeros(Int, K - length(v)))
nonzero(v) = [c > 0 ? c : NaN for c in v]    # no bar for a zero count on the log axis
steps_fig = Plots.bar(
    collect(0:(K - 1)) .- 0.2,
    nonzero(pad(os));
    bar_width = 0.4,
    fillto = 0.8,
    label = "optimistic (some successors)",
    color = :orange,
    yscale = :log10,
    xlabel = "steps to the target",
    ylabel = "cells",
    title = "Cells by number of steps to the target ($tlabel)",
    legend = :topright,
    size = (800, 450),
)
Plots.bar!(
    steps_fig,
    collect(0:(K - 1)) .+ 0.2,
    nonzero(pad(ds));
    bar_width = 0.4,
    fillto = 0.8,
    label = "worst case (the domain)",
    color = :steelblue,
)
Plots.plot!(
    steps_fig;
    ylims = (0.8, 3 * max(maximum(os), maximum(ds))),
    xticks = 0:(K > 25 ? 5 : 1):(K - 1),
)
Plots.savefig(steps_fig, out("_steps.png"))

# ------------------------------------------------------------
# Closed loops from the deepest cells of the domain
# ------------------------------------------------------------
for (run, r) in d["runs"]
    isempty(r.inputs) && (println("run $run: no input"); continue)
    traj = ST.Trajectory(r.states; inputs = r.inputs)
    files = RP.plot_footstep(
        traj,
        disc.state_grid,
        carving.cells;
        x_lb,
        x_ub,
        tstep = disc.tstep,
        geometry,
        obstacle,
        foothold,
        file_prefix = out("_" * run),
        ankle_label = "p₁ (N m s)",
        dynamics,
        subframes = 10,
        title = "Closed loop from the $(replace(run, "_" => " ")) domain cell ($(r.steps) worst-case steps)",
    )
    println(
        "run $run: $(length(r.inputs)) steps, foot from ",
        round.(foot(r.states[1]); digits = 3),
        " to ",
        round.(foot(r.states[end]); digits = 3),
        "; ",
        files,
    )
end

# ------------------------------------------------------------
# Closed loops from the start cells in the domain, side by side
# ------------------------------------------------------------
start_runs = filter(r -> !isempty(r.inputs), get(d, "start_runs", []))
if length(start_runs) > 12
    # A spread of them, from the fewest to the most worst-case steps.
    println("start runs: ", length(start_runs), ", 12 of them drawn")
    sort!(start_runs; by = r -> r.steps)
    start_runs = start_runs[unique(round.(Int, range(1, length(start_runs); length = 12)))]
end
if !isempty(start_runs)
    n = length(start_runs)
    colours = n == 1 ? [:steelblue] : Plots.palette(:viridis, n + 1)[1:n]   # no pale yellow end
    fa = Plots.plot(;
        xlabel = "x (m)",
        ylabel = "z (m)",
        title = "(a) swing foot",
        legend = :bottomleft,
        legendfontsize = 7,
    )
    draw_scene!(fa; xlims = (-0.3, 0.3), ylims = (-0.17, 0.12))
    fq = Plots.plot(;
        xlabel = "time (s)",
        ylabel = "q₁ (rad)",
        title = "(b) ankle angle",
        legend = false,
    )
    fp = Plots.plot(;
        xlabel = "time (s)",
        ylabel = "p₁ (N m s)",
        title = "(c) ankle momentum",
        legend = false,
    )
    for (r, c) in zip(start_runs, colours)
        fs = foot.(r.states)
        ts = (0:(length(r.states) - 1)) .* disc.tstep
        Plots.plot!(
            fa,
            first.(fs),
            last.(fs);
            color = c,
            marker = :circle,
            ms = 2,
            label = @sprintf("p₁ = %.3f, foot x = %.2f", r.start[1], foot(r.start)[1]),
        )
        Plots.plot!(fq, ts, [s[2] for s in r.states]; color = c, marker = :circle, ms = 2)
        Plots.plot!(fp, ts, [s[1] for s in r.states]; color = c, marker = :circle, ms = 2)
        println(
            @sprintf(
                "start run from p₁ = %.4f (%d worst-case steps): %d steps, foot from ",
                r.start[1],
                r.steps,
                length(r.inputs),
            ),
            round.(foot(r.states[1]); digits = 3),
            " to ",
            round.(foot(r.states[end]); digits = 3),
        )
    end
    Plots.savefig(
        Plots.plot(
            fa,
            fq,
            fp;
            layout = Plots.@layout([a{0.45w} [b; c]]),
            size = (1200, 520),
            plot_title = "Closed loops from the start configuration ($tlabel)",
            margin = 4Plots.mm,
        ),
        out("_starts.png"),
    )
    println("saved ", out("_starts.png"))
end

# ------------------------------------------------------------
# The optimistic path from the start, against the RK4 map
# ------------------------------------------------------------
if path !== nothing
    n = length(path.states)
    # Does the RK4 map, with the inputs of the optimistic path, stay in the domain and end in the target?
    X_kept = UT.set_minus(X_box, carving.cells)
    target = RP.target_set(
        disc.state_grid,
        foothold,
        geometry,
        X_box;
        tolerance = d["foothold_tolerance"],
    )
    leaves = findfirst(x -> !(x in X_kept), replay)
    replay_note =
        (
            leaves === nothing ? "stays in the domain" :
            "leaves the domain at step $(leaves - 1)"
        ) *
        (replay[end] in target ? ", ends in the target" : ", does not end in the target")
    println("RK4 replay of the optimistic path: ", replay_note)
    function frame(k)
        left = Plots.plot(; title = "optimistic path: cell $(k - 1) of $(n - 1)")
        draw_scene!(left)
        for j in 1:(k - 1)
            RP.draw_robot!(left, path.states[j]; geometry, alpha = 0.15)
        end
        RP.draw_robot!(left, path.states[k]; geometry)
        right = Plots.plot(; title = "RK4 map, same inputs: step $(k - 1)")
        draw_scene!(right)
        for j in 1:(k - 1)
            RP.draw_robot!(right, replay[j]; geometry, alpha = 0.15)
        end
        RP.draw_robot!(right, replay[k]; geometry)
        u = k < n ? path.inputs[k] : nothing
        sub =
            u === nothing ? "the cells reach the target; the RK4 map $replay_note" :
            "next input (v₂, v₃, v₄) = $(Tuple(round.(u; digits = 1))) rad/s"
        return Plots.plot(
            left,
            right;
            layout = (1, 2),
            size = (1200, 520),
            plot_title = sub,
            plot_titlefontsize = 13,
        )
    end
    anim = Plots.@animate for k in 1:n
        frame(k)
    end
    Plots.gif(anim, out("_optimistic.gif"); fps = 1)
    steps = 0:(n - 1)
    p1 = Plots.plot(
        steps,
        [s[2] for s in path.states];
        marker = :circle,
        color = :orange,
        label = "cells",
        ylabel = "q₁ (rad)",
        title = "ankle along the optimistic path",
    )
    Plots.plot!(
        p1,
        steps,
        [s[2] for s in replay];
        marker = :diamond,
        color = :purple,
        ls = :dash,
        label = "RK4, same inputs",
    )
    p2 = Plots.plot(
        steps,
        [s[1] for s in path.states];
        marker = :circle,
        color = :orange,
        label = "cells",
        ylabel = "p₁ (N m s)",
        xlabel = "step",
    )
    Plots.plot!(
        p2,
        steps,
        [s[1] for s in replay];
        marker = :diamond,
        color = :purple,
        ls = :dash,
        label = "RK4, same inputs",
    )
    Plots.savefig(
        Plots.plot(
            frame(n),
            Plots.plot(p1, p2; layout = (2, 1));
            layout = Plots.@layout([a{0.7w} b]),
            size = (1700, 560),
        ),
        out("_optimistic.png"),
    )
    println(
        "optimistic path: ",
        n - 1,
        " steps; ankle (q₁, p₁) along the cells ",
        [round.((s[2], s[1]); digits = 3) for s in path.states],
        "; along the RK4 map ",
        [round.((s[2], s[1]); digits = 3) for s in replay],
    )
end
println("figures in ", OUTDIR)
