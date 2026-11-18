# Sizing a low rectangle ("guard") on the ground between the footholds of a step, so that the swing
# foot cannot drag along or under the ground, without carving out the start and end cells.
#
# The ground constraint of `RP.infeasible_cells` keeps the cells across the ground by default
# (`ground_ambiguous = :keep`): removing them would remove the footholds too. A kept cell may then put
# the swing foot below the ground, down to the `z_min` of its foot box, so a step could drag the foot
# on or under the floor. A rectangle obstacle `[x_l, x_r] × [z_b, z_t]` between the footholds removes
# every cell whose foot box meets it: the ambiguous cells under it (their box has `z_max ≥ 0`, so it
# meets any rectangle with `z_b ≤ 0 ≤ z_t` that overlaps it along `x`), and the cells where the foot
# is lower than `z_t` above it. The kept cells over the span `[x_l, x_r]` then have the whole foot
# above `z_t` (certified, with a certified `foot_bound`).
#
# For the state box and the resolutions below, the script prints:
#
# * how deep the kept cells put the foot below the ground: the smallest `z_min` over the kept cells,
#   overall, over the step region `[-0.2, 0.2]`, and per 5 cm of `x` (the profile); the dragging that
#   is left outside the span of the rectangle is at most this deep;
# * how far the foot box of a foothold cell reaches towards the middle of the step: for the start
#   foothold `p_s`, the largest `x_max - p_s,x` over its cells; for the end foothold `p_e`, the largest
#   `p_e,x - x_min`. A rectangle with `x_l` above `p_s,x` plus that reach (and `x_r` below `p_e,x`
#   minus it) carves out none of them. The foothold cells are those of `RP.target_cells` (the target
#   set of `RP.step_problem`, sampled) and, as a guaranteed superset, every cell whose foot box
#   contains the foothold;
# * a check: `RP.infeasible_cells` with the largest such rectangle (rounded inwards to the
#   millimetre) and a few heights, against the foothold cells of `RP.target_cells`; also with the
#   rectangle of half-width 0.15 m for comparison.
#
# Set `STATE_BOUNDS`, `RESOLUTIONS`, `FOOTHOLDS = (p_s, p_e)`, `FOOT_BOUND` and `GUARD_HEIGHTS` before
# including this file to change them; `PROFILE_FILE` is the plot of the profile.
#
#   julia --project=problems/BipedRobot tests/step_guard_sizing.jl

using Printf
using StaticArrays
import Dionysos
import LazySets
import Plots

const MP = Dionysos.Mapping

include(joinpath(@__DIR__, "..", "robot_problem.jl"))
import .UnderactuatedRobotProblem as RP

# The knee box of size 1.5: the knees bend backwards, `q₂ ≥ 0` and `q₄ ≤ 0`.
q_lb, q_ub =
    @isdefined(STATE_BOUNDS) ? STATE_BOUNDS :
    (SVector(-1.5, 0.0, -1.5, -1.5), SVector(1.5, 1.5, 1.5, 0.0))
resolutions = @isdefined(RESOLUTIONS) ? RESOLUTIONS : (0.1, 0.05)
p_start, p_end = @isdefined(FOOTHOLDS) ? FOOTHOLDS : (SVector(-0.2, 0.0), SVector(0.2, 0.0))
foot_bound = @isdefined(FOOT_BOUND) ? FOOT_BOUND : :second_order_cell
guard_heights = @isdefined(GUARD_HEIGHTS) ? GUARD_HEIGHTS : (0.02, 0.03, 0.05)
profile_file =
    @isdefined(PROFILE_FILE) ? PROFILE_FILE :
    joinpath(@__DIR__, "ground_cells_plots", "step_guard_profile.png")
const BINS = -0.6:0.05:0.6   # the edges of the bins of the profile along `x`

# One pass over the cells of the box: the kept cells (`z_max ≥ 0`), how deep they put the foot, and
# the reach of the cells whose foot box contains a foothold.
function scan(q_lb, q_ub, h, geometry, foot_bound, p_start, p_end)
    q_grid = MP.GridFree(zero(q_lb), SVector(h, h, h, h))
    Q = LazySets.Hyperrectangle(; low = q_lb, high = q_ub)
    r = SVector(h, h, h, h) ./ 2
    R2, R3 = RP._foot_remainders(geometry, r)
    ncells = nkept = nambiguous = 0
    deepest = (; z = Inf, x = (NaN, NaN), q = zero(q_lb))
    deepest_step = Inf                        # over the kept cells that overlap [p_s,x, p_e,x]
    profile = fill(Inf, length(BINS) - 1)     # per bin of `x`
    reach = [-Inf, -Inf]                      # start, end
    ncontain = [0, 0]
    for pos in MP.get_pos_from_set(q_grid, Q, MP.INNER)
        ncells += 1
        q = SVector{4, Float64}(MP.get_coord_by_pos(q_grid, pos))
        x, z = RP._foot_box(q, r, geometry, foot_bound, R2, R3)
        for (k, p) in enumerate((p_start, p_end))
            if x[1] <= p[1] <= x[2] && z[1] <= p[2] <= z[2]
                ncontain[k] += 1
                reach[k] = max(reach[k], k == 1 ? x[2] - p[1] : p[1] - x[1])
            end
        end
        z[2] < 0 && continue                  # below the ground: removed
        nkept += 1
        z[1] < 0 || continue
        nambiguous += 1
        z[1] < deepest.z && (deepest = (; z = z[1], x, q))
        x[1] <= p_end[1] && p_start[1] <= x[2] && (deepest_step = min(deepest_step, z[1]))
        for b in eachindex(profile)
            x[1] <= BINS[b + 1] && BINS[b] <= x[2] && (profile[b] = min(profile[b], z[1]))
        end
    end
    return (;
        ncells,
        nkept,
        nambiguous,
        deepest,
        deepest_step,
        profile,
        reach,
        ncontain,
        r,
        R3,
    )
end

# The foothold cells of `RP.target_cells` (the target set of the step problem) that lie entirely in the
# box, as 4-D positions, and the reach of their foot boxes towards the middle of the step.
function target_reach(state_grid, X, p, geometry, foot_bound, inward)
    q_grid = MP.GridFree(zero(q_lb), SVector{4}(MP.get_h(state_grid)[2:5]))
    Q = LazySets.Hyperrectangle(; low = q_lb, high = q_ub)
    inner = Set(MP.get_pos_from_set(q_grid, Q, MP.INNER))
    cells = unique(
        pos[2:5] for pos in RP.target_cells(state_grid, p, geometry, X; nsub = 6) if
        pos[2:5] in inner
    )
    r = MP.get_h(q_grid) ./ 2
    best = -Inf
    for pos in cells
        x = RP.swing_foot_range(MP.get_coord_by_pos(q_grid, pos), r, geometry; foot_bound).x
        best = max(best, inward > 0 ? x[2] - p[1] : p[1] - x[1])
    end
    return cells, best
end

# The cells of `cells` (4-D positions) that `carving` removes.
lost(carving, cells) = count(pos -> (0, pos...) in carving.cells, cells)

cm(v) = @sprintf("%.2f cm", 100 * v)
geometry = RP.default_geometry()
profiles = Dict{Float64, Vector{Float64}}()
guards = Dict{Float64, Tuple{Float64, Float64}}()
for h in resolutions
    println(
        "\n==== cell size $h rad, box $(Tuple(q_lb)) – $(Tuple(q_ub)), foot_bound = ",
        repr(foot_bound),
    )
    s = @timed scan(q_lb, q_ub, h, geometry, foot_bound, p_start, p_end)
    res = s.value
    profiles[h] = res.profile
    println(
        "$(res.ncells) cells, $(res.nkept) kept by the ground constraint, $(res.nambiguous) of them ",
        "across the ground (",
        round(s.time; digits = 1),
        " s)",
    )
    println(
        "deepest kept cell: z_min = ",
        cm(res.deepest.z),
        " (foot box x ∈ [",
        cm(res.deepest.x[1]),
        ", ",
        cm(res.deepest.x[2]),
        "], q = ",
        round.(Tuple(res.deepest.q); digits = 2),
        "); over the step region x ∈ [$(p_start[1]), $(p_end[1])]: ",
        cm(res.deepest_step),
        "  (true depth within R₃ = ",
        @sprintf("%.2f mm", 1e3 * res.R3),
        ")",
    )
    print("profile of z_min (kept cells), per 5 cm of x:")
    for b in eachindex(res.profile)
        -0.35 <= BINS[b] < 0.35 || continue
        @printf("\n  x ∈ [%+.2f, %+.2f]: %s", BINS[b], BINS[b + 1], cm(res.profile[b]))
    end
    println()

    # The 5-D grid of the abstraction, with one cell along the ankle axis (the carving ignores it).
    state_grid = MP.GridFree(SVector(0.0, 0, 0, 0, 0), SVector(0.1, h, h, h, h))
    X = LazySets.Hyperrectangle(;
        low = SVector(-0.05, q_lb...),
        high = SVector(0.05, q_ub...),
    )
    start_cells, start_reach =
        target_reach(state_grid, X, p_start, geometry, foot_bound, +1)
    end_cells, end_reach = target_reach(state_grid, X, p_end, geometry, foot_bound, -1)
    println(
        "start foothold $(Tuple(p_start)): $(length(start_cells)) target cells, their foot boxes ",
        "reach up to x = ",
        cm(p_start[1] + start_reach),
        " (+",
        cm(start_reach),
        "); every cell whose box contains it ($(res.ncontain[1])): up to x = ",
        cm(p_start[1] + res.reach[1]),
    )
    println(
        "end foothold $(Tuple(p_end)): $(length(end_cells)) target cells, their foot boxes reach ",
        "down to x = ",
        cm(p_end[1] - end_reach),
        " (-",
        cm(end_reach),
        "); every cell whose box contains it ($(res.ncontain[2])): down to x = ",
        cm(p_end[1] - res.reach[2]),
    )

    # The widest guard that keeps every target cell, rounded inwards to the millimetre, against the
    # target cells; and the guard of half-width 0.15 m.
    x_l = ceil((p_start[1] + start_reach) * 1000 + 1e-9) / 1000
    x_r = floor((p_end[1] - end_reach) * 1000 - 1e-9) / 1000
    guards[h] = (x_l, x_r)
    z_b = floor(res.deepest.z * 1000) / 1000
    for (xl, xr) in ((x_l, x_r), (-0.15, 0.15)), z_t in guard_heights
        guard = LazySets.Hyperrectangle(; low = SVector(xl, z_b), high = SVector(xr, z_t))
        carving = RP.infeasible_cells(
            X,
            state_grid,
            [guard],
            geometry;
            obstacle_name = "guard",
            foot_bound,
        )
        @printf(
            "guard x ∈ [%+.3f, %+.3f], z ∈ [%.3f, %.3f]: removes %d cells (with the ground: %d of %d); target cells lost: start %d / %d, end %d / %d\n",
            xl,
            xr,
            z_b,
            z_t,
            length(last(carving.causes[2])),
            length(carving.cells.positions),
            res.ncells,
            lost(carving, start_cells),
            length(start_cells),
            lost(carving, end_cells),
            length(end_cells),
        )
    end
end

# The profile of the depth of the kept cells along `x`, for every resolution, with the footholds and
# the widest guards.
fig = Plots.plot(;
    xlabel = "swing foot x [m]",
    ylabel = "lowest z of the kept cells [cm]",
    title = "How deep the kept cells put the swing foot (ground_ambiguous = :keep)",
    titlefontsize = 10,
    legend = :outerbottom,
    legend_columns = 3,
    size = (900, 560),
    xlims = (-0.4, 0.4),
)
for (k, h) in enumerate(sort(collect(keys(profiles)); rev = true))
    color = (:steelblue, :darkorange, :seagreen)[mod1(k, 3)]
    Plots.plot!(
        fig,
        collect(BINS),
        100 .* [profiles[h]; profiles[h][end]];
        seriestype = :steppost,
        label = "cell $h rad",
        color,
        linewidth = 2,
    )
    x_l, x_r = guards[h]
    Plots.vspan!(fig, [x_l, x_r]; color, alpha = 0.12, label = "widest guard, cell $h rad")
end
Plots.vline!(
    fig,
    [p_start[1], p_end[1]];
    color = :black,
    linestyle = :dash,
    label = "footholds",
)
Plots.hline!(fig, [0.0]; color = :black, label = "")
mkpath(dirname(profile_file))
Plots.savefig(fig, profile_file)
println("\nProfile saved to ", profile_file)
