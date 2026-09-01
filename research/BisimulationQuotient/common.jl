# Shared setup for the PCLF bisimulation-quotient experiments (HSCC 2027).
#
# Every script in this folder `include`s this file. It is not library code: it holds only the
# boilerplate the experiments genuinely share — package imports, the module aliases, the two
# benchmark problems, the MOI call sequences for building a quotient and synthesising a co-safe
# LTL controller on it, and the figures the paper uses. Anything specific to one experiment
# (the path-complete graph, the specification, the tolerances) stays in that experiment.
#
# `include` is textual, so the `using` lines below also bring the packages into the calling
# script — a script needs no imports of its own beyond `Spot`, which only the LTL experiments
# pay for.

using StaticArrays
using LinearAlgebra
using JuMP
using Clarabel
using JLD2
using CDDLib

import HybridSystems
using LazySets
using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

const PCQ = AB.PCLFBisimulationQuotient
const PCLF = UT.PathCompleteFramework

# ---------------------------------------------------------
# Persisting a solved optimizer
# ---------------------------------------------------------
# Constructing a quotient costs minutes, so the expensive experiments save theirs and reload it
# while iterating on the figures. `*.jld2` is gitignored.

function export_optimizer_jld2(optimizer, filename::AbstractString)
    jldopen(filename, "w") do file
        file["format_version"] = 1
        return file["optimizer"] = optimizer
    end
    return nothing
end

function import_optimizer_jld2(filename::AbstractString)
    return jldopen(filename, "r") do file
        v = file["format_version"]
        v == 1 || error("Unsupported optimizer file format_version=$v")
        return file["optimizer"]
    end
end

# ---------------------------------------------------------
# The two benchmark problems
# ---------------------------------------------------------

"""
Two contracting modes on a square, with two rectangular observation regions.

Used by the De Bruijn and custom-graph experiments. Returns the switched system, the
`BisimulationQuotientProblem`, and the regions, which the LTL specifications refer to by name.
"""
function two_mode_problem()
    A1 = @SMatrix [0.70 0.10; 0.00 0.65]
    A2 = @SMatrix [0.60 -0.15; 0.10 0.55]
    f = HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)])

    X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
    R1 = LazySets.Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5])
    R2 = LazySets.Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5])

    problem = PR.BisimulationQuotientProblem(f, X, [R1, R2])
    return (; f, problem, X, R1, R2)
end

"""
The two-mode system of the paper's observer-graph study, on a square of half-width `p`.

`rotation` turns the state-space polytope by that angle, which is what distinguishes the
rotated-state-space experiment from the PCLF/CLF comparison.
"""
function observer_graph_problem(; p::Float64 = 1.7, rotation::Float64 = 0.0)
    A1 = (1.0 / 10.0) * [1.5519 0.4474; 7.6412 7.4716]
    A2 = (1.0 / 10.0) * [0.4750 9.1755; 1.8955 0.1850]
    f = HybridSystems.discreteswitchedsystem([A1, A2])

    Rot = [cos(rotation) -sin(rotation); sin(rotation) cos(rotation)]
    normals = [Rot * n for n in ([1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0])]
    X = LazySets.HPolytope([LazySets.HalfSpace(n, p) for n in normals])

    R1 = LazySets.HPolytope([
        LazySets.HalfSpace([1.0, 0.0], 1.5),
        LazySets.HalfSpace([-1.0, 0.0], -0.8),
        LazySets.HalfSpace([0.0, 1.0], 1.5),
        LazySets.HalfSpace([0.0, -1.0], -0.8),
    ])

    problem = PR.BisimulationQuotientProblem(f, X, [R1])
    return (; f, problem, X, R1)
end

"""
Example 3.1 of Gol, Ding, Lazar & Belta, *Finite Bisimulations for Switched Linear Systems*
(arXiv:1208.5471), with the three observation regions used throughout this folder.

Shared so that our certificate and theirs are compared on the same problem: the working set and
the regions must match, or the two quotients partition different sets and their cell counts are
not comparable.

Their published certificate is `V(x) = ‖Lx‖_∞` with `L = gol_lazar_belta_L()`, contracting at
`ρ = 0.94`.
"""
function gol_lazar_belta_problem(; p::Float64 = 5.9)
    A1 = @SMatrix [-0.65 0.32; -0.42 -0.92]
    A2 = @SMatrix [0.65 0.32; -0.42 -0.92]
    f = HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)])

    X = LazySets.HPolytope([
        LazySets.HalfSpace([1.0, 0.0], p),
        LazySets.HalfSpace([-1.0, 0.0], p),
        LazySets.HalfSpace([0.0, 1.0], p),
        LazySets.HalfSpace([0.0, -1.0], p),
    ])

    R1 = LazySets.HPolytope(
        [-0.9869 -0.1615; -0.0931 0.9957; 0.9659 0.2587; 0.0825 -0.9966],
        [6.6767, 9.2315, 2.3700, -5.9038],
    )
    R2 = LazySets.HPolytope(
        [0.9993 0.0363; -0.7743 -0.6329; 0.5463 0.8376],
        [-2.1809, 6.3754, -4.8983],
    )
    R3 = LazySets.HPolytope(
        [-0.9946 -0.1041; 0.5277 0.8494; 0.9999 0.0146; -0.1191 -0.9929],
        [-5.5771, 5.3510, 9.1600, 6.2406],
    )

    problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])
    return (; f, problem, X, R1, R2, R3)
end

"""
The polytopic Lyapunov function published with the example above: `V(x) = ‖Lx‖_∞`, `ρ = 0.94`.
"""
gol_lazar_belta_L() = [
    -0.0625 1.0
    0.6815 1.0
    0.9947 0.6868
    0.9947 -0.0678
]

"""
The four-node observer graph and the conic partition both experiments on it share.
"""
function observer_graph_pclf(f; max_iter::Int = 1000)
    graph = PCLF.edgeList_to_LabDigraph([
        (1, 2, 1),
        (2, 1, 1),
        (2, 4, 1),
        (2, 3, 1),
        (3, 4, 1),
        (4, 3, 2),
        (4, 4, 2),
        (4, 1, 2),
    ])

    rays = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-1.0, 1.0], [-1.0, 0.0]]
    cones = [hcat(rays[i], rays[i + 1]) for i in 1:(length(rays) - 1)]
    partitions = Dict(v => cones for v in 1:4)

    sdp_optimizer =
        JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => max_iter)
    return PCLF.compute_polyhedral_pieces_pclf(
        f,
        graph,
        sdp_optimizer,
        partitions;
        MLF = true,
    )
end

# `Gmats` rotates the template of each node, which is what lets a path-complete family beat a
# common Lyapunov function on this system.
function rotation_templates(nodes; θ::Float64 = π / 6, mode::Symbol = :rotation, n::Int = 2)
    Rot = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    if mode == :rotation
        return Dict(v => Rot for v in nodes)
    elseif mode == :identity
        return Dict(v => Matrix{Float64}(I, n, n) for v in nodes)
    elseif mode == :alternating
        return Dict(v => (iseven(sum(v)) ? Rot : Matrix{Float64}(I, n, n)) for v in nodes)
    else
        error(
            "Unknown template mode = $mode; expected :rotation, :identity or :alternating",
        )
    end
end

# ---------------------------------------------------------
# The two solver call sequences
# ---------------------------------------------------------

"""
Build the bisimulation quotient of `problem` under `pclf`.

Returns the optimizer alongside the quotient and its terminal set `D`, because the optimizer is
what `export_optimizer_jld2` persists.
"""
function build_quotient(
    problem,
    pclf;
    atol::Float64 = 1e-4,
    level_tol::Union{Nothing, Float64} = nothing,
    max_levels::Union{Nothing, Int} = nothing,
    max_slices::Int = 10,
    print_level::Int = 1,
    backend = nothing,
)
    optimizer = MOI.instantiate(PCQ.OptimizerBisimulationQuotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), atol)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), max_slices)
    isnothing(level_tol) ||
        MOI.set(optimizer, MOI.RawOptimizerAttribute("level_tol"), level_tol)
    isnothing(max_levels) ||
        MOI.set(optimizer, MOI.RawOptimizerAttribute("max_levels"), max_levels)
    isnothing(backend) ||
        MOI.set(optimizer, MOI.RawOptimizerAttribute("polyhedra_backend"), backend)

    MOI.optimize!(optimizer)
    println(
        "Construction time = ",
        MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec")),
    )

    quotient = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
    D = MOI.get(optimizer, MOI.RawOptimizerAttribute("D"))
    PCQ.print_bisimulation_stats(quotient)
    return (; optimizer, quotient, D)
end

"""
Synthesise a co-safe LTL controller on `quotient` and simulate it from `x0`.

`regions` maps each atomic proposition to its set and `ap_to_obs` to the observation label the
quotient carries for it; the latter is what synthesis actually reads. `early_stop` decides what
`success` means -- see the `OptimizerCoSafeLTLOnQuotient` docstring.
"""
function synthesize_cosafe_ltl(
    f,
    quotient,
    spec,
    regions::Dict{Symbol},
    ap_to_obs::Dict{Symbol, Int},
    x0;
    early_stop::Bool = false,
    print_level::Int = 1,
    N::Int = 50,
)
    initial_set = LazySets.Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]])
    problem = PR.CoSafeLTLProblem(
        f,
        initial_set,
        spec,
        regions,
        Dict{Symbol, Any}(ap => MP.INNER for ap in keys(regions)),
    )

    optimizer = MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"), quotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("ap_to_obs"), ap_to_obs)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), early_stop)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)
    MOI.optimize!(optimizer)

    controller = PCQ.solve_concrete_problem(optimizer)
    controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
    uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))

    mem0 = PCQ.initial_controller_memory(optimizer, x0)
    sim = PCQ.simulate_closed_loop(f, controller, x0, mem0; N = N)

    return (;
        optimizer,
        controller,
        controllable_set,
        uncontrollable_set,
        X = sim.X,
        U = sim.U,
        M = sim.M,
    )
end

# ---------------------------------------------------------
# Figures
# ---------------------------------------------------------

# Node keys are `Int` for a hand-written graph and `Tuple` for a De Bruijn one; both order
# naturally, and sorting by `string` would put node 10 before node 2.
quotient_nodes(quotient) = sort(collect(keys(quotient.slices)))

"""
One panel per automaton node, showing every quotient cell of that node.
"""
function plot_cells_by_node(quotient; nodes = quotient_nodes(quotient))
    fig = plot(; layout = (1, length(nodes)), aspect_ratio = :equal, legend = false)
    for (i, nd) in enumerate(nodes)
        plot!(fig[i], quotient; what = :states, node = nd, show_contours = false)
        title!(fig[i], "Node $nd")
    end
    return fig
end

"""
The winning region alone, with the closed-loop trajectory over it.
"""
function plot_controllable_set(
    quotient,
    controllable_set,
    X_seq,
    title;
    xlims = (-4.5, 4.5),
    ylims = (-4.5, 4.5),
)
    fig = plot(; aspect_ratio = :equal, legend = false, title = title)
    plot!(
        fig,
        quotient;
        what = :states,
        state_ids = controllable_set,
        show_contours = false,
        user_color = :green,
        fillalpha = 1.0,
    )
    plot!(fig, ST.Trajectory(X_seq); label = "Trajectory")
    xlims!(fig, xlims...)
    ylims!(fig, ylims...)
    return fig
end

# Green over red rather than red over green: the winning region is the result, so it goes on top.
function _plot_winning_losing!(
    panel,
    quotient,
    controllable_set,
    uncontrollable_set,
    X_seq;
    kw...,
)
    plot!(
        panel,
        quotient;
        what = :states,
        state_ids = uncontrollable_set,
        show_contours = false,
        user_color = :red,
        fillalpha = 1.0,
        kw...,
    )
    plot!(
        panel,
        quotient;
        what = :states,
        state_ids = controllable_set,
        show_contours = false,
        user_color = :green,
        fillalpha = 1.0,
        kw...,
    )
    plot!(panel, ST.Trajectory(X_seq); label = "Trajectory")
    return panel
end

"""
Winning and losing cells per node, then all nodes together, printing each node's volume.
"""
function plot_partition_by_node(
    quotient,
    controllable_set,
    uncontrollable_set,
    X_seq;
    nodes = quotient_nodes(quotient),
    xlims = (-4.5, 4.5),
    ylims = (-4.5, 4.5),
    backend = CDDLib.Library(),
)
    L = length(nodes)
    fig = plot(; layout = (1, L + 1), aspect_ratio = :equal, legend = false)

    for (i, nd) in enumerate(nodes)
        winning_here = PCQ.state_ids_in_node(quotient, nd; state_ids = controllable_set)
        volume = PCQ.get_volume(quotient, winning_here; backend = backend)
        println("Volume of controllable set in node $nd = ", volume)

        title!(fig[i], "Node $nd")
        xlims!(fig[i], xlims...)
        ylims!(fig[i], ylims...)
        _plot_winning_losing!(
            fig[i],
            quotient,
            controllable_set,
            uncontrollable_set,
            X_seq;
            node = nd,
        )
    end

    println(
        "Volume of controllable set over all nodes = ",
        PCQ.get_volume(quotient, controllable_set; backend = backend),
    )
    title!(fig[L + 1], "All states")
    xlims!(fig[L + 1], xlims...)
    ylims!(fig[L + 1], ylims...)
    _plot_winning_losing!(fig[L + 1], quotient, controllable_set, uncontrollable_set, X_seq)

    return fig
end
