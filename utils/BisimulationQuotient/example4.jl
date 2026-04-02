using StaticArrays
using LinearAlgebra
using JuMP
using Clarabel

import HybridSystems
using LazySets
using Plots
using Spot
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

const PCLF = UT.PathCompleteFramework
# =========================================================
# Generic experiment for De Bruijn graphs
# =========================================================

function plot_partition_by_node(
    bisimulation,
    controllable_set,
    uncontrollable_set,
    X_seq;
    xlims = (-4.5, 4.5),
    ylims = (-4.5, 4.5),
)
    nodes = sort(collect(keys(bisimulation.slices)))
    L = length(nodes)

    fig = plot(; layout = (1, L + 1), aspect_ratio = :equal, legend = false)

    for (i, nd) in enumerate(nodes)
        ctrl_node = AB.PCLFBisimulationQuotient.state_ids_in_node(
            bisimulation,
            nd;
            state_ids = controllable_set,
        )
        Vctrl_node = AB.PCLFBisimulationQuotient.get_volume(bisimulation, ctrl_node)
        println("Volume of controllable set in Node $nd = ", Vctrl_node)

        title!(fig[i], "Node $nd")
        xlims!(fig[i], xlims...)
        ylims!(fig[i], ylims...)

        plot!(
            fig[i],
            bisimulation;
            what = :states,
            state_ids = controllable_set,
            node = nd,
            show_contours = false,
            user_color = :green,
            fillalpha = 1.0,
        )
        plot!(
            fig[i],
            bisimulation;
            what = :states,
            state_ids = uncontrollable_set,
            node = nd,
            show_contours = false,
            user_color = :red,
            fillalpha = 1.0,
        )
        plot!(fig[i], ST.Trajectory(X_seq); label = "Trajectory")
    end

    Vctrl = AB.PCLFBisimulationQuotient.get_volume(bisimulation, controllable_set)
    println("Volume of controllable set in all states = ", Vctrl)

    title!(fig[L + 1], "All states")
    xlims!(fig[L + 1], xlims...)
    ylims!(fig[L + 1], ylims...)

    plot!(
        fig[L + 1],
        bisimulation;
        what = :states,
        state_ids = uncontrollable_set,
        show_contours = false,
        user_color = :red,
        fillalpha = 1.0,
    )
    plot!(
        fig[L + 1],
        bisimulation;
        what = :states,
        state_ids = controllable_set,
        show_contours = false,
        user_color = :green,
        fillalpha = 1.0,
    )
    plot!(fig[L + 1], ST.Trajectory(X_seq); label = "Trajectory")

    return fig
end

function plot_controllable_set(
    bisimulation,
    controllable_set,
    X_seq,
    φ_str;
    xlims = (-4.5, 4.5),
    ylims = (-4.5, 4.5),
)
    fig = plot(; aspect_ratio = :equal, legend = false)
    plot!(fig; title = φ_str)
    plot!(
        fig,
        bisimulation;
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

function make_Gmats(nodes; θ = π / 6, mode = :rotation, n = 2)
    R = [
        cos(θ) -sin(θ)
        sin(θ) cos(θ)
    ]

    if mode == :rotation
        return Dict(v => R for v in nodes)
    elseif mode == :identity
        return Dict(v => Matrix{Float64}(I, n, n) for v in nodes)
    elseif mode == :alternating
        return Dict(v => (iseven(sum(v)) ? R : Matrix{Float64}(I, n, n)) for v in nodes)
    else
        error("Unknown Gmats mode = $mode")
    end
end

function run_example_debruijn(;
    k::Int = 1,
    dual::Bool = true,
    gmode::Symbol = :rotation,
    θ::Float64 = π / 6,
    max_iter::Int = 1000,
    max_levels::Int = 100,
    max_slices::Int = 8,
    atol::Float64 = 1e-3,
    verbose_bisim::Bool = true,
    early_stop::Bool = false,
    print_level::Int = 1,
)
    # ---------------------------------------------------------
    # Define a stable switched system
    # ---------------------------------------------------------
    A1 = @SMatrix [
        0.70 0.10
        0.00 0.65
    ]

    A2 = @SMatrix [
        0.60 -0.15
        0.10 0.55
    ]

    f = HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)])

    # ---------------------------------------------------------
    # Define problem
    # ---------------------------------------------------------
    X = Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
    R1 = Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5])
    R2 = Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5])
    observation_regions = [R1, R2]
    problem = PR.BisimulationQuotientProblem(f, X, observation_regions)

    # ---------------------------------------------------------
    # Build PCLF on De Bruijn graph
    # ---------------------------------------------------------
    G = PCLF.generate_DeBruijn_edges(2, k; dual = dual)
    nodes = sort(collect(G.verts))

    sdp_optimizer =
        JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => max_iter)

    Gmats = make_Gmats(nodes; θ = θ, mode = gmode, n = 2)

    pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
        f,
        G,
        sdp_optimizer;
        Gmats = Gmats,
        MLF = true,
        verbose = false,
    )

    println()
    println("------------------------------------------")
    println("k = $k")
    println("Graph nodes = ", nodes)
    println("Number of graph nodes = ", length(nodes))
    println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

    # ---------------------------------------------------------
    # Compute bisimulation quotient
    # ---------------------------------------------------------
    optimizer_bisim =
        MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)

    MOI.set(
        optimizer_bisim,
        MOI.RawOptimizerAttribute("bisimulation_quotient_problem"),
        problem,
    )
    MOI.set(optimizer_bisim, MOI.RawOptimizerAttribute("pclf"), pclf)
    MOI.set(optimizer_bisim, MOI.RawOptimizerAttribute("verbose"), verbose_bisim)
    MOI.set(optimizer_bisim, MOI.RawOptimizerAttribute("atol"), atol)
    MOI.set(optimizer_bisim, MOI.RawOptimizerAttribute("max_levels"), max_levels)
    MOI.set(optimizer_bisim, MOI.RawOptimizerAttribute("max_slices"), max_slices)

    MOI.optimize!(optimizer_bisim)

    construction_time =
        MOI.get(optimizer_bisim, MOI.RawOptimizerAttribute("construction_time_sec"))
    println("Construction time = ", construction_time)

    bisimulation =
        MOI.get(optimizer_bisim, MOI.RawOptimizerAttribute("bisimulation_quotient"))
    D = MOI.get(optimizer_bisim, MOI.RawOptimizerAttribute("D"))

    AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

    # ---------------------------------------------------------
    # Plot all quotient cells by node
    # ---------------------------------------------------------
    fig_cells = plot(; layout = (1, length(nodes)), aspect_ratio = :equal, legend = false)

    for (i, nd) in enumerate(nodes)
        plot!(fig_cells[i], bisimulation; what = :states, node = nd, show_contours = false)
        title!(fig_cells[i], "Node $nd")
    end
    display(fig_cells)

    # ---------------------------------------------------------
    # CoSafe LTL control synthesis on the quotient
    # ---------------------------------------------------------
    φ = ltl"F(R1 & F(D))"
    φ_str = string(φ)

    x0 = SVector(2.3, 1.5)
    _I_ = Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]])

    prob = PR.CoSafeLTLProblem(
        f,
        _I_,
        φ,
        Dict(:D => D, :R1 => R1, :R2 => R2),
        Dict{Symbol, Any}(:D => MP.INNER, :R1 => MP.INNER, :R2 => MP.INNER),
        true,
    )

    optimizer_ltl =
        MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerCoSafeLTLOnQuotient)
    MOI.set(optimizer_ltl, MOI.RawOptimizerAttribute("concrete_problem"), prob)
    MOI.set(optimizer_ltl, MOI.RawOptimizerAttribute("bisimulation_quotient"), bisimulation)
    MOI.set(
        optimizer_ltl,
        MOI.RawOptimizerAttribute("ap_to_obs"),
        Dict(:D => -1, :R1 => 1, :R2 => 2),
    )
    MOI.set(optimizer_ltl, MOI.RawOptimizerAttribute("early_stop"), early_stop)
    MOI.set(optimizer_ltl, MOI.RawOptimizerAttribute("print_level"), print_level)

    MOI.optimize!(optimizer_ltl)

    concrete_controller = AB.PCLFBisimulationQuotient.solve_concrete_problem(optimizer_ltl)
    controllable_set = MOI.get(optimizer_ltl, MOI.RawOptimizerAttribute("controllable_set"))
    uncontrollable_set =
        MOI.get(optimizer_ltl, MOI.RawOptimizerAttribute("uncontrollable_set"))

    mem0 = AB.PCLFBisimulationQuotient.initial_controller_memory(optimizer_ltl, x0)

    sim = AB.PCLFBisimulationQuotient.simulate_closed_loop(
        f,
        concrete_controller,
        x0,
        mem0;
        N = 50,
    )

    X_seq = sim.X
    U_seq = sim.U
    M_seq = sim.M

    println("X_seq = ", X_seq)
    println("U_seq = ", U_seq)
    println("M_seq = ", M_seq)

    # ---------------------------------------------------------
    # Plot controllable set only
    # ---------------------------------------------------------
    fig_ctrl = plot_controllable_set(
        bisimulation,
        controllable_set,
        X_seq,
        φ_str;
        xlims = (-4.5, 4.5),
        ylims = (-4.5, 4.5),
    )
    display(fig_ctrl)

    # ---------------------------------------------------------
    # Plot by node + all states
    # ---------------------------------------------------------
    fig_partition = plot_partition_by_node(
        bisimulation,
        controllable_set,
        uncontrollable_set,
        X_seq;
        xlims = (-4.5, 4.5),
        ylims = (-4.5, 4.5),
    )
    return display(fig_partition)
end

# =========================================================
# Run examples
# For M = 2 modes, De Bruijn graphs have:
# k = 0 -> 1 node
# k = 1 -> 2 nodes
# k = 2 -> 4 nodes
# =========================================================

run_example_debruijn(; k = 0, gmode = :rotation)
run_example_debruijn(; k = 1, gmode = :rotation)
run_example_debruijn(; k = 2, gmode = :rotation)
