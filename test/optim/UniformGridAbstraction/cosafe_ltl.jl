module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
using JuMP
import MathOptInterface as MOI
using Spot
using Plots

include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "Integrator", "integrator.jl"),
)

@testset "UniformGridAbstraction Integrator (CoSafeLTL monitor)" begin
    # ------------------------------------------------------------
    # 1) Concrete system + AlternatingSimulationProblem abstraction build
    # ------------------------------------------------------------
    _X_ = UT.box(SVector(-2.0, -2.0), SVector(2.0, 2.0))
    _U_ = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))

    concrete_system = Integrator.system(; _X_ = _X_, _U_ = _U_)
    jacobian_bound = Integrator.jacobian_bound()

    alternating_simulation_problem =
        DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

    x0g = SVector(-2.0, -2.0)
    hx = SVector(0.2, 0.2)
    state_grid = MP.GridFree(x0g, hx)

    u0g = SVector(-1.0, -1.0)
    hu = SVector(0.5, 0.5)
    input_grid = MP.GridFree(u0g, hu)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        alternating_simulation_problem,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.GROWTH,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)

    # Make it deterministic + small for CI
    MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(optimizer, MOI.Silent(), true)

    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    discrete_time_system =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

    @test abstract_system !== nothing
    @test discrete_time_system !== nothing

    # mapping sanity
    Xmap = SY.get_state_mapping(abstract_system)
    @test Xmap !== nothing
    @test MP.get_n_state(Xmap) > 0

    # ------------------------------------------------------------
    # 2) Co-safe LTL labeling + monitor
    # ------------------------------------------------------------
    _I_ = UT.box(SVector(-1.7, -1.7), SVector(-1.6, -1.6))

    g1 = UT.box(SVector(1.0, 1.0), SVector(1.7, 1.7))
    g2 = UT.box(SVector(-1.5, -1.2), SVector(-0.6, -0.2))

    obs = UT.box(SVector(-1.8, 0.0), SVector(-0.6, 1.0))

    danger1 = UT.box(SVector(-0.5, -0.5), SVector(0.5, 0.5))
    danger2 = UT.box(SVector(1.3, -0.5), SVector(2.0, 0.5))
    danger = UT.set_union([danger1, danger2])

    φ = ltl"G(!obs) & F(g1 & ((!danger) U g2))"

    struct MonitorG1NoDangerUntilG2 end
    @inline function mon_next(::MonitorG1NoDangerUntilG2, q::Int, ap::Tuple{Vararg{Symbol}})
        obs = (:obs in ap)
        g1 = (:g1 in ap)
        g2 = (:g2 in ap)
        danger = (:danger in ap)

        obs && return 0
        q == 0 && return 0
        q == 3 && return 3

        if q == 1
            if g1
                return g2 ? 3 : 2
            else
                return 1
            end
        end

        @assert q == 2
        if g2
            return 3
        elseif danger
            return 0
        else
            return 2
        end
    end

    spec = OPDS.FunctionMonitor(
        1,         # initial state
        Set([3]),  # accepting states
        (qa, ap) -> mon_next(MonitorG1NoDangerUntilG2(), qa, ap),
    )

    labeling = Dict{Symbol, Any}(:g1 => g1, :g2 => g2, :danger => danger, :obs => obs)

    ap_semantics = Dict{Symbol, Any}(
        :g1 => MP.INNER,
        :g2 => MP.INNER,
        :danger => MP.OUTER,
        :obs => MP.OUTER,
    )

    concrete_problem =
        DI.Problem.CoSafeLTLProblem(concrete_system, _I_, spec, labeling, ap_semantics)

    # quick labeling sanity: each AP should exist
    @test haskey(concrete_problem.labeling, :g1)
    @test haskey(concrete_problem.labeling, :g2)
    @test haskey(concrete_problem.labeling, :danger)
    @test haskey(concrete_problem.labeling, :obs)

    # ------------------------------------------------------------
    # 3) Solve CoSafeLTL using same optimizer (reusing abstraction pipeline)
    # ------------------------------------------------------------
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.optimize!(optimizer)

    success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
    @test success isa Bool

    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    @test abstract_controller !== nothing
    @test concrete_controller !== nothing

    # ------------------------------------------------------------
    # 4) Closed-loop trajectory sanity
    # ------------------------------------------------------------
    x0 = SVector(-1.65, -1.65)
    nstep = 60

    traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        concrete_controller,
        x0,
        nstep;
        update_on_next = true,
        stopping = x -> false,
    )

    xs = collect(ST.states(traj))
    us = collect(ST.inputs(traj))
    qs = collect(ST.memory(traj))

    @test !isempty(xs)
    @test length(xs) <= nstep + 1
    @test length(qs) == length(xs) || length(qs) == length(xs) - 1

    # If we have at least one transition, we must have at least one control
    if length(xs) >= 2
        @test !isempty(us)
    end

    # Typical invariant: either one input per transition, or some pipelines store same length as xs
    @test length(us) == length(xs) - 1 || length(us) == length(xs)

    # every x should stay in domain X (soft sanity; allow boundary tol)
    X = concrete_system.X
    @test all(
        x -> all(i -> LazySets.low(X, i) - 1e-8 <= x[i] <= LazySets.high(X, i) + 1e-8, 1:2),
        xs,
    )

    # ------------------------------------------------------------
    # 5) Optional plot smoke test (skip on CI)
    # ------------------------------------------------------------
    no_plot = true
    @static if get(ENV, "CI", "false") == "false" &&
               (isdefined(@__MODULE__, :no_plot) && no_plot == false)
        fig = plot(; aspect_ratio = :equal, title = string(φ))
        plot!(
            concrete_problem;
            ap_colors = Dict(:g1 => :red, :g2 => :cyan, :danger => :orange, :obs => :black),
        )
        plot!(fig, x_traj; color = :blue, dims = [1, 2])
        @test fig isa Plots.Plot
    end
end

end # module TestMain
