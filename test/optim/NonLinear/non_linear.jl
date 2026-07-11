module TestNonLinearBemporadMorari

include("../GolLazarBelta/solvers.jl")
include("../../../problems/NonLinear/non_linear.jl")

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim

using LinearAlgebra, Test
using MathematicalSystems, HybridSystems
using FillArrays

@testset "NonLinear BemporadMorari" begin
    T = Float64

    # Get the nonlinear system from the problem definition
    nl_prob = NonLinear.problem()
    nl_sys = nl_prob.system

    # Wrap f_eval into a (x, u) -> next_state function (fixing w = 0)
    w_zero = zeros(nl_sys.nw)
    f = (x, u) -> nl_sys.f_eval(x, u, w_zero)

    nlmap = BlackBoxControlDiscreteSystem(f, nl_sys.nx, nl_sys.nu)

    # State domain as a box (automatic conversion to MOI constraints)
    pX = UT.box([-5.0, -5.0], [5.0, 5.0])

    automaton = GraphAutomaton(1)
    add_transition!(automaton, 1, 1, 1)

    sys = HybridSystem(
        automaton,
        [ConstrainedContinuousIdentitySystem(2, pX)],
        Fill(nlmap, 1),
        Fill(ControlledSwitching(), 1),
    )

    # Use Ipopt for NLP (nonlinear constraints)
    nlp_solver =
        optimizer_with_attributes(cont_solver.optimizer_constructor, MOI.Silent() => true)

    @testset "Depth $N" for (N, x_0) in [(1, T[-2.0, -2.0]), (3, T[-3.0, -3.0])]
        q_0 = 1
        q_T = 1

        state_cost = Fill(UT.ZeroFunction(), nmodes(sys))
        transition_cost = UT.QuadraticFunction(Matrix{T}(I, 2, 2))
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
            "continuous_solver" => nlp_solver,
            "mixed_integer_solver" => miqp_solver,
            "indicator" => false,
            "print_level" => 0,
        )

        optimizer = MOI.instantiate(algo)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
        MOI.optimize!(optimizer)

        @test MOI.get(optimizer, MOI.TerminationStatus()) in
              [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED]

        xu = MOI.get(optimizer, ST.ContinuousTrajectoryAttribute())
        @test length(xu.x) == N
        @test length(xu.u) == N

        # Verify dynamics: x_next ≈ f(x_prev, u) using the original f_eval
        for t in 1:N
            x_prev = t == 1 ? x_0 : xu.x[t - 1]
            x_next = xu.x[t]
            u_t = xu.u[t]
            f_val = nl_sys.f_eval(x_prev, u_t, w_zero)
            @test x_next ≈ f_val atol = 1e-4
        end

        obj = MOI.get(optimizer, MOI.ObjectiveValue())
        @test isfinite(obj)
        @test obj >= -1e-6
    end
end

end # module
