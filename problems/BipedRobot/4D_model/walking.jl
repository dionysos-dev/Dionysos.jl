# ============================================================
# Walking: the foot strike as a hybrid transition
# ============================================================
#
# The state is *role-relative* — `(hip_stance, knee_stance, hip_swing, knee_swing)`
# with the stance foot at the origin — so a step is one continuous swing phase
# followed by a discrete strike. Writing the old swing chain from its own foot,
# once that foot becomes the new stance foot, gives the reset exactly:
#
#     (θ₁,θ₂,θ₃,θ₄) ↦ (θ₃,θ₄,θ₁,θ₂)
#
# with no translation term in the state — the frame hop lives in the world, which
# is not a state. Two consequences: the reset is an exact lattice automorphism
# when the four axes share a grid step (so the abstraction stays exact through the
# strike), and both stance phases are *the same system*, so they share one
# abstraction. Which physical leg is planted is the mode, not a state.

"""
    leg_swap(θ)

The strike reset `(θ₁,θ₂,θ₃,θ₄) ↦ (θ₃,θ₄,θ₁,θ₂)`: the swing leg becomes the stance
leg. An involution, and an exact grid automorphism when all four axes share a step.
"""
leg_swap(θ) = SVector(θ[3], θ[4], θ[1], θ[2])

"""
    strike_guard(X, state_grid, geometry; ground_band, min_advance, grounded_left_foot = true)

Cells where the strike is available: the swing foot is within `ground_band` of the
ground and at least `min_advance` ahead of the stance foot.

A band, not the surface `foot_y = 0`: an equality guard has measure zero, so no
cell lies *inside* it and the guard would discretize to nothing. `ground_band`
should be at least the Lipschitz deviation of a cell
([`swing_foot_deviation_bound`](@ref)), which is the resolution at which "the foot
is on the ground" is decidable at all.

`min_advance` is the stride the guard demands. Synthesis minimises the number of
steps, so the gait strides exactly as little as the guard allows — it is a
specification, not an emergent property.
"""
function strike_guard(
    X::LazySets.AbstractHyperrectangle,
    state_grid::MP.GridFree,
    geometry::Geometry;
    ground_band::Real = swing_foot_deviation_bound(geometry, MP.get_h(state_grid) ./ 2),
    min_advance::Real = 0.05,
    grounded_left_foot::Bool = true,
)
    return MP.cells_where(state_grid, X) do pos
        foot = swing_foot_position(
            SVector{4}(MP.get_coord_by_pos(state_grid, pos)),
            geometry,
            grounded_left_foot,
        )
        return abs(foot[2]) <= ground_band && foot[1] >= min_advance
    end
end

"""
    walking_hybrid_system(concrete_system, guard)

Two-mode hybrid system for walking: mode 1 = left foot planted, mode 2 = right
foot planted, both carrying `concrete_system` — *the same object*, so the two
modes share one abstraction — and both transitions guarded by `guard` with the
[`leg_swap`](@ref) reset.

`concrete_system` must be the **continuous-time** form (see [`system`](@ref)):
the hybrid machinery integrates a mode's dynamics itself over the mode's time
step, so it needs the vector field, not the one-step map.

The switch is offered to the synthesis rather than forced: the impossibility of
sinking through the ground is already carried by the carved domain, and in
quasi-static walking *when* to put the foot down is genuinely a control decision.
"""
function walking_hybrid_system(concrete_system, guard)
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    HybridSystems.add_transition!(automaton, 2, 1, 2)
    strike = ST.GuardedResetMap(guard, leg_swap)
    return HybridSystems.HybridSystem(
        automaton,
        [concrete_system, concrete_system],
        [strike, strike],
        [HybridSystems.AutonomousSwitching(), HybridSystems.AutonomousSwitching()],
    )
end

"""
    walking_problem(; geometry, disc, x_bar, x0, min_advance, obstacle)

Everything a walking run needs, built once: the carved continuous-time system,
the strike guard, the two-mode hybrid system, and the reach-the-guard
specification whose solution the [`WalkingController`](@ref) chains into a gait.

Returns `(; hs, guard, system, domain, X, x0, problem, optimizer_kwargs)`.
"""
function walking_problem(;
    geometry::Geometry = default_geometry(),
    disc = default_discretization(; dx = 0.1, tstep = 0.1, speed_levels = 1),
    x_bar::Real = 0.6,
    x0::SVector{4} = SVector(0.2, 0.0, -0.2, 0.0),
    min_advance::Real = 0.15,
    obstacle = nothing,
)
    X = LazySets.Hyperrectangle(;
        low = SVector(-x_bar, -x_bar, -x_bar, -x_bar),
        high = SVector(x_bar, x_bar, x_bar, x_bar),
    )
    domain = domain_for(disc; x_bar = x_bar)
    removed = infeasible_cells(X, disc.state_grid, obstacle, geometry)
    concrete_system = system(;
        tstep = disc.tstep,
        domain = domain,
        state_grid = disc.state_grid,
        removed_cells = removed,
        continuous_time = true,
    )
    guard = strike_guard(X, disc.state_grid, geometry; min_advance = min_advance)
    hs = walking_hybrid_system(concrete_system, guard)

    # Reach the guard, in either mode: the gait is the recurrence of that reach
    # (see `gait_recurrence`), not a safety property — standing still would
    # satisfy safety, since `u = 0` is an input.
    target = PR.HybridSpec(Dict(1 => PR.StateSpec(guard), 2 => PR.StateSpec(guard)))
    problem = PR.OptimalControlProblem(
        hs,
        (x0, 1),
        target,
        nothing,
        (aug, u) -> 1.0,
        PR.Infinity(),
    )

    optimizer_kwargs = Dict(
        "state_grid" => disc.state_grid,
        "input_grid" => disc.input_grid,
        "time_step" => disc.tstep,
        "approx_mode" => DI.Optim.Abstraction.UniformGridAbstraction.CENTER_SIMULATION,
        "print_level" => 0,
    )

    return (; hs, guard, system = concrete_system, domain, X, x0, problem, optimizer_kwargs)
end

"""
    solve_walking(w; print_level = 0)

Build the abstraction of `w = walking_problem(...)` and synthesize the
reach-the-guard controller, returning
`(; optimizer, abstract_system, controllable, success, walker)` where `walker`
is the [`WalkingController`](@ref) that chains footsteps.

Both modes carry the same system, so the abstraction underneath is built once.
"""
function solve_walking(w; print_level::Int = 0)
    optimizer = MOI.instantiate(DI.Optim.Abstraction.HybridSystemAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), w.problem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("optimizer_list"),
        Function[
            () -> MOI.instantiate(DI.Optim.Abstraction.UniformGridAbstraction.Optimizer),
            () -> MOI.instantiate(DI.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        ],
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
        [copy(w.optimizer_kwargs), copy(w.optimizer_kwargs)],
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    walker = WalkingController(
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller")),
        w.guard,
        w.optimizer_kwargs["state_grid"],
        w.hs,
        abstract_system,
    )
    return (;
        optimizer,
        abstract_system,
        controllable = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set")),
        success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success")),
        walker,
    )
end

"""
    walk(w, solution, nstep)

Closed-loop walking trajectory: `(aug_states, inputs)` with the augmented state
`(x, mode)`, the strikes appearing as switching labels among the inputs.
"""
walk(w, solution, nstep::Integer) =
    DI.Optim.Abstraction.HybridSystemAbstraction.get_closed_loop_trajectory(
        w.hs,
        solution.walker,
        [w.optimizer_kwargs["time_step"], w.optimizer_kwargs["time_step"]],
        (w.x0, 1),
        nstep,
    )

"""
    gait_recurrence(abstract_system, guard, state_grid, controllable) -> (checked, failed)

The walking certificate. With `G` the guard and `Reach(G)` the states from which
`G` is reachable, `π(G) ⊆ Reach(G)` means: from any post-strike state, reach the
guard again, strike, repeat — so the robot can strike **forever**, established by
one finite reachability synthesis rather than a Büchi condition.

Returns how many post-strike states were checked and how many are *not*
recurrent; `failed == 0` is the certificate. A non-zero count names a genuine
gap — the remedy is a wider guard or a finer grid, not a change of model.
"""
function gait_recurrence(abstract_system, guard, state_grid::MP.GridFree, controllable)
    checked, failed = 0, 0
    for pos in guard, next_mode in (2, 1)
        swapped = leg_swap(SVector{4}(MP.get_coord_by_pos(state_grid, pos)))
        q = SY.get_abstract_state(abstract_system, (swapped, next_mode))
        (q === nothing || q <= 0) && continue
        checked += 1
        q in controllable || (failed += 1)
    end
    return checked, failed
end

"""
    WalkingController(reach_controller, guard, state_grid, hs, abstract_system)

Chains footsteps: strike whenever the guard is reached, and otherwise follow
`reach_controller`, the controller synthesized to *reach* the guard.

Nothing more is needed to walk forever. The reachability controller is defined on
every state from which the guard is reachable, and the recurrence certificate
([`gait_recurrence`](@ref)) is exactly the statement that the post-strike image
lands back in that set — so the same controller applies after every strike, with
no re-synthesis. This is the 4-D counterpart of the 6D model's two-step walking
controller: a phase wrapper around already-synthesized pieces.

The strike labels are read from the abstraction's input map rather than
formatted here: the spelling of a switching input is a contract owned by the
abstraction and parsed back by the simulation.
"""
struct WalkingController{C, G, S} <: ST.AbstractContinuousController
    reach_controller::C
    guard::G
    state_grid::S
    strike_label::Dict{Int, String} # stance mode → label of the strike leaving it
end

function WalkingController(reach_controller, guard, state_grid, hs, abstract_system)
    input_map = SY.get_global_input_map(abstract_system)
    labels = Dict{Int, String}()
    for (transition_id, transition) in
        enumerate(collect(HybridSystems.transitions(hs.automaton)))
        source = HybridSystems.source(hs.automaton, transition)
        labels[source] = SY.get_switch_label(input_map, transition_id)
    end
    return WalkingController(reach_controller, guard, state_grid, labels)
end

ST.controller_kind(::WalkingController) = ST.StaticKind()

# `aug_state` is the hybrid `(x, mode)`.
_in_guard(ctrl::WalkingController, aug_state) =
    MP.get_pos_by_coord(ctrl.state_grid, SVector{4}(aug_state[1])) in ctrl.guard

function ST.is_defined(ctrl::WalkingController, mem, aug_state)
    _in_guard(ctrl, aug_state) && haskey(ctrl.strike_label, aug_state[end]) && return true
    return ST.is_defined(ctrl.reach_controller, mem, aug_state)
end

function ST.output_control(ctrl::WalkingController, mem, aug_state)
    if _in_guard(ctrl, aug_state)
        label = get(ctrl.strike_label, aug_state[end], nothing)
        label === nothing || return label
    end
    return ST.output_control(ctrl.reach_controller, mem, aug_state)
end

"""
    physical_pose(θ, left_stance) -> (θ_drawn, grounded_left_foot)

Arguments for the kinematics from a role-relative state and the physical stance
side.

The state always carries the stance leg in `θ[1:2]`, whereas
[`cartesian_coordinates`](@ref) reads the stance from `θ[3:4]` when the *right*
foot is planted. Passing the state through unchanged would therefore plant the
swinging leg every other step — the posture comes out mirrored, and the walk
looks like it steps backwards half the time. Swapping the pairs realigns the two
conventions.
"""
physical_pose(θ::SVector{4}, left_stance::Bool) =
    left_stance ? (θ, true) : (leg_swap(θ), false)

"""
    walk_world_offsets(aug_states, inputs, geometry)

Cumulative world position of the stance foot along a walking trajectory: the
frame hops forward by the landing position of the swing foot at every strike.
The offset is *driver-side bookkeeping*, not a state — the model is
stance-relative, which is exactly why one abstraction serves every step.

Returns one offset per state, and the physical `grounded_left_foot` flag per
state (it alternates at each strike, so the legs keep their identity on screen).
"""
function walk_world_offsets(aug_states, inputs, geometry::Geometry = default_geometry())
    offsets = zeros(length(aug_states))
    left_stance = trues(length(aug_states))
    offset = 0.0
    left = true
    for k in eachindex(aug_states)
        offsets[k] = offset
        left_stance[k] = left
        k <= length(inputs) || continue
        if inputs[k] isa AbstractString
            # The strike: the new stance foot is where the swing foot just landed.
            foot = swing_foot_position(SVector{4}(aug_states[k][1]), geometry, true)
            offset += foot[1]
            left = !left
        end
    end
    return offsets, left_stance
end
