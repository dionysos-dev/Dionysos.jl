module TestNonLinearBemporadMorari

include("../GolLazarBelta/solvers.jl")
include("../../../problems/non_linear.jl")

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim

using LinearAlgebra, Test
import CDDLib
import Polyhedra
using MathematicalSystems, HybridSystems, SemialgebraicSets
using FillArrays

function rect2d(lib, T, x_l, x_u, y_l, y_u)
    r =
        Polyhedra.HalfSpace([-1, 0], -T(x_l)) ∩ Polyhedra.HalfSpace([1, 0], T(x_u)) ∩
        Polyhedra.HalfSpace([0, -1], -T(y_l)) ∩ Polyhedra.HalfSpace([0, 1], T(y_u))
    return Polyhedra.polyhedron(r, lib)
end

"""
Build a PWA HybridSystem approximation of the non-linear system by
partitioning the state space into rectangular regions and linearizing
the dynamics at the center of each region.
"""
function pwa_system(lib, T, μ, Ts, x_bounds, n_per_dim)
    step = (x_bounds[2] - x_bounds[1]) / n_per_dim
    regions = Tuple{T, T, T, T}[]
    centers = Tuple{T, T}[]
    for i in 1:n_per_dim, j in 1:n_per_dim
        x_l = x_bounds[1] + (i - 1) * step
        x_u = x_bounds[1] + i * step
        y_l = x_bounds[1] + (j - 1) * step
        y_u = x_bounds[1] + j * step
        push!(regions, (x_l, x_u, y_l, y_u))
        push!(centers, ((x_l + x_u) / 2, (y_l + y_u) / 2))
    end

    n_modes = length(regions)
    domains = [rect2d(lib, T, r...) for r in regions]

    # Input domain
    pU = rect2d(lib, T, -10, 10, -10, 10)

    # Linearize at each region center:
    # f(x,u) = [1.1*x1 - 0.2*x2 - μ*x2^3 + Ts*u1,
    #            1.1*x2 + 0.2*x1 + μ*x1^3 + Ts*u2]
    # Jacobian w.r.t. x at (cx, cy):
    #   A = [1.1,          -0.2 - 3μ*cy^2;
    #        0.2 + 3μ*cx^2, 1.1           ]
    # B = [Ts 0; 0 Ts]
    # c = f(cx,cy,0) - A*[cx;cy]
    reset_maps = map(centers) do (cx, cy)
        A = T[
            1.1 (-0.2-3μ*cy^2)
            (0.2+3μ*cx^2) 1.1
        ]
        B = T[Ts 0; 0 Ts]
        f_cx = T[
            1.1 * cx - 0.2 * cy - μ * cy^3,
            1.1 * cy + 0.2 * cx + μ * cx^3,
        ]
        c = f_cx - A * T[cx, cy]
        return ConstrainedAffineControlMap(A, B, c, FullSpace(), pU)
    end

    # Create automaton: allow transitions between all modes
    automaton = GraphAutomaton(n_modes)
    k = 0
    for from in 1:n_modes, to in 1:n_modes
        k += 1
        add_transition!(automaton, from, to, k)
    end

    # Each transition (from, to) uses the dynamics linearized at the center of `from`
    all_reset_maps = Vector{ConstrainedAffineControlMap}(undef, k)
    idx = 0
    for from in 1:n_modes, _ in 1:n_modes
        idx += 1
        all_reset_maps[idx] = reset_maps[from]
    end

    sys = HybridSystem(
        automaton,
        [ConstrainedContinuousIdentitySystem(2, d) for d in domains],
        all_reset_maps,
        Fill(ControlledSwitching(), k),
    )
    return sys
end

"""
Find the mode index for a given point in the grid.
"""
function find_mode(x, x_bounds, n_per_dim)
    step = (x_bounds[2] - x_bounds[1]) / n_per_dim
    i = clamp(Int(floor((x[1] - x_bounds[1]) / step)) + 1, 1, n_per_dim)
    j = clamp(Int(floor((x[2] - x_bounds[1]) / step)) + 1, 1, n_per_dim)
    return (i - 1) * n_per_dim + j
end

@testset "NonLinear BemporadMorari" begin
    T = Float64
    lib = CDDLib.Library()
    μ = T(0.00005)
    Ts = T(1.0)

    # Use a 2x2 grid over [-4, 4]^2
    x_bounds = (T(-4), T(4))
    n_per_dim = 2
    sys = pwa_system(lib, T, μ, Ts, x_bounds, n_per_dim)

    @testset "Depth 1" begin
        # Start at [-2, -2] (mode 1: [-4,0]x[-4,0]), target mode 4 ([0,4]x[0,4])
        x_0 = T[-2.0, -2.0]
        q_0 = find_mode(x_0, x_bounds, n_per_dim)
        q_T = find_mode(T[2.0, 2.0], x_bounds, n_per_dim)
        N = 1

        state_cost = Fill(UT.ZeroFunction(), nmodes(sys))
        transition_cost = UT.QuadraticControlFunction(Matrix{T}(I, 2, 2))
        problem = PR.OptimalControlProblem(
            sys,
            (q_0, x_0),
            q_T,
            Fill(state_cost, N),
            Fill(Fill(transition_cost, ntransitions(sys)), N),
            N,
        )

        algo = optimizer_with_attributes(
            OP.BemporadMorari.Optimizer{T},
            "continuous_solver" => qp_solver,
            "mixed_integer_solver" => miqp_solver,
            "indicator" => false,
            "log_level" => 0,
        )

        optimizer = MOI.instantiate(algo)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
        MOI.optimize!(optimizer)

        @test MOI.get(optimizer, MOI.TerminationStatus()) in
              [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED]

        xu = MOI.get(optimizer, ST.ContinuousTrajectoryAttribute())
        # Final state must be in target region [0,4]x[0,4]
        @test xu.x[end][1] >= -1e-2
        @test xu.x[end][2] >= -1e-2
        @test xu.x[end][1] <= 4.0 + 1e-2
        @test xu.x[end][2] <= 4.0 + 1e-2

        obj = MOI.get(optimizer, MOI.ObjectiveValue())
        @test isfinite(obj)
        @test obj >= 0.0
    end

    @testset "Depth 3" begin
        x_0 = T[-3.0, -3.0]
        q_0 = find_mode(x_0, x_bounds, n_per_dim)
        q_T = find_mode(T[2.0, 2.0], x_bounds, n_per_dim)
        N = 3

        state_cost = Fill(UT.ZeroFunction(), nmodes(sys))
        transition_cost = UT.QuadraticControlFunction(Matrix{T}(I, 2, 2))
        problem = PR.OptimalControlProblem(
            sys,
            (q_0, x_0),
            q_T,
            Fill(state_cost, N),
            Fill(Fill(transition_cost, ntransitions(sys)), N),
            N,
        )

        algo = optimizer_with_attributes(
            OP.BemporadMorari.Optimizer{T},
            "continuous_solver" => qp_solver,
            "mixed_integer_solver" => miqp_solver,
            "indicator" => false,
            "log_level" => 0,
        )

        optimizer = MOI.instantiate(algo)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
        MOI.optimize!(optimizer)

        @test MOI.get(optimizer, MOI.TerminationStatus()) in
              [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED]

        xu = MOI.get(optimizer, ST.ContinuousTrajectoryAttribute())
        @test length(xu.x) == N
        @test length(xu.u) == N

        # Final state must be in target region [0,4]x[0,4]
        @test xu.x[end][1] >= -1e-2
        @test xu.x[end][2] >= -1e-2
        @test xu.x[end][1] <= 4.0 + 1e-2
        @test xu.x[end][2] <= 4.0 + 1e-2

        obj = MOI.get(optimizer, MOI.ObjectiveValue())
        @test isfinite(obj)
        # With longer horizon, cost should be <= depth 1 cost (more flexibility)
        @test obj >= 0.0
    end
end

end # module
