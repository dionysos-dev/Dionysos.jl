# # 4-D velocity-controlled biped: certified multi-step walking
#
# The single certified footstep of `biped_4d_velocity.jl` becomes *walking* by
# modelling the foot strike as a hybrid transition: one continuous mode (the leg
# swing) and a discrete switch whose guard is "the swing foot is on the ground,
# ahead" and whose reset swaps the legs.
#
# Three facts make this small rather than expensive:
#
# 1. **The state does not grow.** It stays role-relative —
#    `(hip_stance, knee_stance, hip_swing, knee_swing)` with the stance foot at the
#    origin — so the horizontal position never becomes a state. It is recovered
#    driver-side by accumulating the step length (`RP.walk_world_offsets`); the
#    velocity `ẋ` has no place at all in a quasi-static model.
# 2. **One abstraction serves both stance phases.** Both modes carry the very
#    same system, so `HybridSystemAbstraction` builds the abstraction once and
#    reuses it — the leg swap lives entirely in the reset map.
# 3. **The reset is exact.** `(θ₁,θ₂,θ₃,θ₄) ↦ (θ₃,θ₄,θ₁,θ₂)` maps grid cells onto
#    grid cells when the four axes share a step, so the abstraction stays exact
#    through the strike, not just during the swing.
#
# The specification is a **recurrence**, not a safety property: staying inside the
# domain forever is trivially satisfied by standing still (`u = 0` is an input).
# What must hold is that the strike can be taken again and again — and a finite
# computation certifies that infinite behaviour. With `G` the guard and `Reach(G)`
# the states from which `G` is reachable, `π(G) ⊆ Reach(G)` means: from any
# post-strike state, reach the guard again, strike, repeat. One reachability
# synthesis, no Büchi machinery.

using StaticArrays
using Dionysos
import MathOptInterface as MOI
import HybridSystems
import LazySets
import Plots

const DI = Dionysos
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "BipedRobot",
        "4D_model",
        "robot_problem.jl",
    ),
)
import .RobotProblem as RP

# ------------------------------------------------------------
# 1) The walking model
# ------------------------------------------------------------

geometry = RP.default_geometry()
disc = RP.default_discretization(; dx = 0.1, tstep = 0.1, speed_levels = 1)

x_bar = 0.6
X = LazySets.Hyperrectangle(;
    low = SVector(-x_bar, -x_bar, -x_bar, -x_bar),
    high = SVector(x_bar, x_bar, x_bar, x_bar),
)
domain = RP.RobotDomainConfig(;
    x_lb = SVector{4}(LazySets.low(X)),
    x_ub = SVector{4}(LazySets.high(X)),
    u_lb = SVector(-disc.u_max, -disc.u_max, -disc.u_max, -disc.u_max),
    u_ub = SVector(disc.u_max, disc.u_max, disc.u_max, disc.u_max),
)

println("Carving the domain (ground constraint only — flat terrain)…")
removed = RP.infeasible_cells(X, disc.state_grid, nothing, geometry)

# A mode carries the vector field `ẋ = u`: the hybrid machinery integrates the
# dynamics itself over the mode's time step, so it needs the field, not the
# one-step map. Both forms give the same abstraction — integrating a constant
# field is exact — but only this one simulates correctly.
concrete_system = RP.system(;
    tstep = disc.tstep,
    domain = domain,
    state_grid = disc.state_grid,
    removed_cells = removed,
    continuous_time = true,
)

# `min_advance` is the stride the guard demands. The synthesis minimises the
# number of steps, so it strides exactly as little as the guard allows: asking
# for 5 cm gives a certified but mincing gait. 15 cm is a stride worth watching,
# and still well inside the leg's ~37 cm reach.
guard = RP.strike_guard(X, disc.state_grid, geometry; min_advance = 0.15)
hs = RP.walking_hybrid_system(concrete_system, guard)
println("guard cells: ", length(guard))

# ------------------------------------------------------------
# 2) Synthesis: reach the guard
# ------------------------------------------------------------

x0 = SVector(0.2, 0.0, -0.2, 0.0)
target = PR.HybridSpec(Dict(1 => PR.StateSpec(guard), 2 => PR.StateSpec(guard)))
problem =
    PR.OptimalControlProblem(hs, (x0, 1), target, nothing, (aug, u) -> 1.0, PR.Infinity())

kwargs = Dict(
    "state_grid" => disc.state_grid,
    "input_grid" => disc.input_grid,
    "time_step" => disc.tstep,
    "approx_mode" => AB.UniformGridAbstraction.CENTER_SIMULATION,
    "print_level" => 0,
)

optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("optimizer_list"),
    Function[
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
    ],
)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
    [copy(kwargs), copy(kwargs)],
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

@time MOI.optimize!(optimizer)
println(
    "reach-the-guard success: ",
    MOI.get(optimizer, MOI.RawOptimizerAttribute("success")),
)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
controllable = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))

# ------------------------------------------------------------
# 3) The recurrence certificate
# ------------------------------------------------------------

n_post_strike, n_not_recurrent = 0, 0
for pos in guard, (mode, next_mode) in ((1, 2), (2, 1))
    swapped = RP.leg_swap(SVector{4}(MP.get_coord_by_pos(disc.state_grid, pos)))
    q = SY.get_abstract_state(abstract_system, (swapped, next_mode))
    (q === nothing || q <= 0) && continue
    global n_post_strike += 1
    q in controllable || (global n_not_recurrent += 1)
end
println(
    "recurrence: $(n_post_strike - n_not_recurrent)/$n_post_strike post-strike states ",
    "can reach the guard again",
)
n_not_recurrent == 0 && println("  ⟹ the robot can strike forever")

# ------------------------------------------------------------
# 4) Walk
# ------------------------------------------------------------

walker = RP.WalkingController(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller")),
    guard,
    disc.state_grid,
)

aug_xs, us = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
    hs,
    walker,
    [disc.tstep, disc.tstep],
    (x0, 1),
    60,
)

strikes = count(u -> u isa AbstractString, us)
offsets, left_stance = RP.walk_world_offsets(aug_xs, us, geometry)
println("steps simulated: ", length(aug_xs) - 1, " with ", strikes, " strikes")
println("distance walked: ", round(offsets[end]; digits = 3), " m")

# ------------------------------------------------------------
# 5) Visualization
# ------------------------------------------------------------
#
# The state is stance-relative, so the world position is *this* accumulator, not
# something the abstraction ever represented.

xlims = (-0.4, offsets[end] + 0.4)
ylims = (-0.05, 0.5)

# Silhouettes along the walk: the robot advances across the ground.
fig = Plots.plot(;
    aspect_ratio = :equal,
    legend = false,
    xlims = xlims,
    ylims = ylims,
    xlabel = "x (m)",
    ylabel = "y (m)",
    title = "Certified walking — $strikes footsteps",
)
Plots.plot!(fig, [xlims[1], xlims[2]], [0.0, 0.0]; color = :black, lw = 1)
# One silhouette per strike (plus the start): the stride sequence, not a blur.
strike_frames = vcat(1, [k for k in eachindex(us) if us[k] isa AbstractString])
for (i, k) in enumerate(strike_frames)
    pose, glf = RP.physical_pose(SVector{4}(aug_xs[k][1]), left_stance[k])
    RP.draw_robot!(
        fig,
        pose;
        geometry = geometry,
        grounded_left_foot = glf,
        origin = offsets[k],
        alpha = 0.3 + 0.7 * (i / length(strike_frames)),
    )
end
# Footprints: where each strike planted the new stance foot.
footprints = unique(offsets)
Plots.scatter!(
    fig,
    footprints,
    zeros(length(footprints));
    marker = :vline,
    markersize = 8,
    color = :green,
)
png_file = joinpath(@__DIR__, "biped_4d_walking.png")
Plots.savefig(fig, png_file)
println("saved ", png_file)

# The animation, drawn frame by frame so each pose sits at its world position.
anim = Plots.@animate for k in 1:length(aug_xs)
    frame = Plots.plot(;
        aspect_ratio = :equal,
        legend = false,
        xlims = xlims,
        ylims = ylims,
        xlabel = "x (m)",
        title = "step $(count(u -> u isa AbstractString, us[1:max(k - 1, 0)])) / $strikes",
    )
    Plots.plot!(frame, [xlims[1], xlims[2]], [0.0, 0.0]; color = :black, lw = 1)
    Plots.scatter!(
        frame,
        footprints,
        zeros(length(footprints));
        marker = :vline,
        markersize = 8,
        color = :green,
    )
    pose, glf = RP.physical_pose(SVector{4}(aug_xs[k][1]), left_stance[k])
    RP.draw_robot!(
        frame,
        pose;
        geometry = geometry,
        grounded_left_foot = glf,
        origin = offsets[k],
    )
    frame
end
gif_file = joinpath(@__DIR__, "biped_4d_walking.gif")
Plots.gif(anim, gif_file; fps = 8)
println("saved ", gif_file)
