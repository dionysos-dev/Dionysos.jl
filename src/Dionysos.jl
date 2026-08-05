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
    add_mode!,
    add_transition!,
    set_role!,
    STATE,
    INPUT,
    CLOCK,
    PARAMETER,
    DISTURBANCE,
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
    add_mode!,
    add_transition!,
    set_role!,
    @mode,
    @specification

# ----- CSV functions for optional dependencies ---------

function export_controller_csv(args...; kwargs...)
    return error("export_controller_csv requires CSV.jl and DataFrames.jl.")
end

function import_controller_csv(args...; kwargs...)
    return error("import_controller_csv requires CSV.jl and DataFrames.jl.")
end

# ----- Spot functions for optional dependencies ---------

function spot_stepper(args...; kwargs...)
    return error(
        "spot_stepper is available only when Spot.jl is installed and loaded with `using Spot`.",
    )
end

# ----- Animation functions for optional dependencies ---------

function animate_trajectory_dashboard end

# ----- RigidBodyDynamics functions for optional dependencies ---------

function animate_mechanism_trajectory end

# ----- Makie functions for optional dependencies ---------

"""
    plot_lifted_bisimulation!(ax, bisimulation; kwargs...)

Draw the lifted quotient states of a PCLF bisimulation quotient as stacked 3D polygons on a
Makie 3D axis `ax`, one horizontal layer per automaton node (its z-height given by `node_z`).

Available only when a Makie backend is loaded (`using GLMakie`, `using CairoMakie`, …); the
method is provided by `DionysosMakieExt`. Keyword arguments select which states to draw
(`state_ids`), how to colour them (`color_by`, `node_colors`), and the outline/opacity style.
"""
function plot_lifted_bisimulation! end

"""
    plot_lifted_trajectory!(ax, bisimulation, state_seq, memory_seq; kwargs...)

Overlay a closed-loop trajectory on the lifted 3D quotient view drawn by
`plot_lifted_bisimulation!`: each planar state in `state_seq` is lifted to the z-height of the
quotient node named by the corresponding entry of `memory_seq`.

Available only when a Makie backend is loaded; the method is provided by `DionysosMakieExt`.
"""
function plot_lifted_trajectory! end

end
