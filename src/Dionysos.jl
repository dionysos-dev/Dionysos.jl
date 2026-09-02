module Dionysos

include("utils/utils.jl")
include("system/system.jl")
include("problem/problem.jl")
include("mapping/mapping.jl")
include("symbolic/symbolic.jl")
include("optim/optim.jl")

# ----- JuMP/MOI front-end ---------
import JuMP
import LazySets

include("wrapper/wrapper.jl")

# Re-export the front-end vocabulary from the top level. `using` (rather than `const` aliases)
# keeps these the *same* bindings as in `Wrapper`, so their docstrings resolve here too.
using .Wrapper:
    Optimizer,
    simulate,
    ∂,
    Δ,
    final,
    start,
    Start,
    Final,
    Always,
    EventuallyAlways,
    Label,
    Guard,
    add_mode!,
    add_transition!,
    set_role!,
    STATE,
    INPUT,
    CLOCK,
    PARAMETER,
    DISTURBANCE,
    CONTINUOUS,
    DISCRETE,
    @mode,
    @specification
# The role constants stay qualified (`Dionysos.STATE`): the bare names are too generic to put
# in a user's namespace.
export ∂,
    Δ,
    final,
    start,
    Start,
    Final,
    Always,
    EventuallyAlways,
    Label,
    Guard,
    add_mode!,
    add_transition!,
    set_role!,
    @mode,
    @specification

# ----- CSV functions for optional dependencies ---------

"""
    export_controller_csv(controller, mapping, filename)

Write a synthesized controller to `filename` as CSV: one row per state of its domain, giving
the cell coordinate and the control input applied there. The controller is read through the
`System.AbstractController` protocol, so any controller type can be exported.

Available only when CSV.jl and DataFrames.jl are loaded; the method is provided by
`DionysosCSVExt`.
"""
function export_controller_csv(args...; kwargs...)
    return error("export_controller_csv requires CSV.jl and DataFrames.jl.")
end

"""
    import_controller_csv(filename, mapping)

Read back a controller written by [`export_controller_csv`](@ref), reattaching it to the
grid `mapping` it was discretized on.

Available only when CSV.jl and DataFrames.jl are loaded; the method is provided by
`DionysosCSVExt`.
"""
function import_controller_csv(args...; kwargs...)
    return error("import_controller_csv requires CSV.jl and DataFrames.jl.")
end

# ----- Spot functions for optional dependencies ---------

"""
    spot_stepper(formula)

Translate a co-safe LTL formula into the deterministic automaton that drives
[`Problem.CoSafeLTLProblem`](@ref) synthesis: a stepper that advances the automaton state as
the closed loop visits labelled regions.

Available only when Spot.jl is installed and loaded (`using Spot`); the method is provided by
`DionysosSpotExt`.
"""
function spot_stepper(args...; kwargs...)
    return error(
        "spot_stepper is available only when Spot.jl is installed and loaded with `using Spot`.",
    )
end

# ----- Animation functions for optional dependencies ---------

"""
    animate_trajectory_dashboard(system_plot!, traj; kwargs...)

Animate a closed-loop trajectory as a multi-panel dashboard: a drawing of the physical
system next to the state, input and — for a hybrid trajectory — mode channels, all advancing
together.

`system_plot!` draws one frame of the plant, and is called as `system_plot!(fig, x, u)`, or
`system_plot!(fig, x, u, mode)` when `traj` carries modes. The `problems/<Name>/` modules
provide one for each benchmark problem. `traj` is a [`System.Trajectory`](@ref) *with
inputs* — the one [`simulate`](@ref Dionysos.Wrapper.simulate) returns.

Returns a `Plots.Animation` when `filename` is omitted, and the filename otherwise. To
display it in a notebook or in the documentation, wrap the animation: `gif(anim; fps = 10)`
— a bare `Plots.Animation` has no `show` method and renders as nothing.

Available only when Plots.jl is loaded (`using Plots`); the method is provided by
`DionysosPlotsExt`.

# Keyword arguments
- `xdims`, `udims`: which coordinates to plot; one dimension gives a time series, two give a
  projection.
- `Δt`: physical time per frame, used for the time axes.
- `fps`: playback speed when writing a file.
- `frame_step`: render every `frame_step`-th frame. This is what shortens render time —
  `fps` only changes playback speed.
- `filename`: write a `.gif` or `.mp4` instead of returning the animation.
- axis labels and limits: `title`, `xlabel_state`, `ylabel_state`, `xlabel_input`,
  `ylabel_input`, `ylabel_mode`, `xlims_system`, `ylims_system`, `xlims_state`,
  `ylims_state`, `xlims_input`, `ylims_input`.
- `show_full_state_traj`, `show_full_input_traj`: draw the whole trajectory faintly behind
  the part played so far.
- `state_background!`: optional `(p_state, x) -> p_state` closure drawn on the state panel
  each frame before the trajectory — e.g. the slice of a carved state region at the current
  state.

```julia
using Plots
anim = Dionysos.animate_trajectory_dashboard(
    system_plot!, trajectory; xdims = (1, 2), udims = (1,), Δt = 0.1, frame_step = 2,
)
gif(anim; fps = 10)
```
"""
function animate_trajectory_dashboard end

# ----- RigidBodyDynamics functions for optional dependencies ---------

"""
    animate_mechanism_trajectory(urdf, traj; joint_names, configuration, Δt, kwargs...)

Replay a closed-loop trajectory on a 3-D mechanism described by a URDF file, using
`RigidBodyDynamics` and `MeshCat`. `configuration` maps a state of `traj` to the joint
configuration vector of the joints named in `joint_names`, and `Δt` is the physical time per
step.

Available only when RigidBodyDynamics.jl, MeshCat.jl and MeshCatMechanisms.jl are loaded;
the method is provided by `DionysosRigidBodyDynamicsExt`.
"""
function animate_mechanism_trajectory end

# ----- Makie functions for optional dependencies ---------

"""
    plot_augmented_bisimulation!(ax, bisimulation; kwargs...)

Draw the node-augmented quotient states of a PCLF bisimulation quotient as stacked 3D polygons on a
Makie 3D axis `ax`, one horizontal layer per automaton node (its z-height given by `node_z`).

Available only when a Makie backend is loaded (`using GLMakie`, `using CairoMakie`, …); the
method is provided by `DionysosMakieExt`. Keyword arguments select which states to draw
(`state_ids`), how to colour them (`color_by`, `node_colors`), and the outline/opacity style.

By default the cells are batched into one mesh per colour and a single stroke for every
outline, since a quotient of a few thousand states would otherwise become as many plot
objects. `merge_plots = false` restores one mesh per polytope, which matters only when the
drawing order of individual overlapping cells does.
"""
function plot_augmented_bisimulation! end

"""
    plot_augmented_trajectory!(ax, bisimulation, state_seq, memory_seq; kwargs...)

Overlay a closed-loop trajectory on the augmented 3D quotient view drawn by
`plot_augmented_bisimulation!`: each planar state in `state_seq` is raised to the z-height of the
quotient node named by the corresponding entry of `memory_seq`.

Available only when a Makie backend is loaded; the method is provided by `DionysosMakieExt`.
"""
function plot_augmented_trajectory! end

end
