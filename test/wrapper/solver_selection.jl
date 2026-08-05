module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Dionysos
using JuMP
using Symbolics
using MathOptSymbolicAD
import MathOptInterface as MOI

# Phase 4: choosing the solver, and the attribute replay that makes the choice possible after
# the options have already been set.

const WR = Dionysos.Wrapper

function reach_model()
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0, start = -0.75)
    @variable(model, -1.0 <= u <= 1.0)
    @constraint(model, ∂(x) == u)
    @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    set_attribute(model, "time_step", 1.0)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(model, "print_level", 0)
    return model
end

@testset "the default solver is inferred, and reported" begin
    model = reach_model()
    @test get_attribute(model, "solver") === nothing      # nothing set ⇒ inferred
    optimize!(model)
    @test backend(model).inner isa AB.UniformGridAbstraction.Optimizer
    @test is_solved_and_feasible(model)
end

@testset "an explicit solver is honoured" begin
    model = reach_model()
    set_attribute(model, "solver", AB.UniformGridAbstraction.Optimizer)
    @test get_attribute(model, "solver") === AB.UniformGridAbstraction.Optimizer
    optimize!(model)
    @test is_solved_and_feasible(model)
end

@testset "options set before the solver are replayed onto it" begin
    # The grids and time step are set first, then the solver is swapped in. Without replay the
    # new solver would come up bare and fail with "Please set the `state_grid`".
    model = reach_model()
    set_attribute(model, "solver", AB.UniformGridAbstraction.Optimizer)

    @test get_attribute(model, "time_step") == 1.0
    @test get_attribute(model, "state_grid") isa MP.Grid

    optimize!(model)
    @test is_solved_and_feasible(model)
end

@testset "unknown options are still rejected at set time" begin
    # Recording options must not turn a typo into a silent no-op discovered much later.
    model = reach_model()
    @test_throws MOI.UnsupportedAttribute set_attribute(model, "no_such_attribute", 1)
end

@testset "an unsupported solver/problem pair is reported by the front-end" begin
    problem = WR.build_problem(backend(reach_model()).ir, (x, u) -> u; time_step = nothing)

    # The uniform grid abstraction handles every problem type the wrapper can build.
    @test WR.supports_problem(AB.UniformGridAbstraction.Optimizer, typeof(problem))
    # An unrelated optimizer declares nothing, so the pair is rejected by name.
    @test !WR.supports_problem(OP.BemporadMorari.Optimizer, typeof(problem))
    err = try
        WR._check_supported(OP.BemporadMorari.Optimizer, problem)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("OptimalControlProblem", err.msg)
end

end # module TestMain
