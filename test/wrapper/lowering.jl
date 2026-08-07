module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Dionysos
using JuMP
using Symbolics
using MathOptSymbolicAD
import MathOptInterface as MOI
import LazySets

# The core lowering contract of the JuMP front-end: how a model becomes a `(system, problem)`
# pair. Started life as the characterisation net for the front-end refactor, so several
# assertions guard historical defects and are tagged FIXED-Ln.
#
# `test/optim/jump_frontend.jl` covers the happy path end to end; this file inspects the
# lowered pair directly, which is both far faster and far more precise.
# `test/wrapper/specifications.jl` covers the specification markers and problem inference.

const WR = Dionysos.Wrapper

# `direct_model` puts our optimizer straight behind JuMP with no caching layer, so the
# constraints land in it as written and the lowered problem can be built without paying
# for an abstraction.
function lowered_problem(build!)
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    build!(model)
    opt = backend(model)
    return WR.lower(opt), opt
end

@testset "lowered problem: shape of the reach-avoid pipeline" begin
    problem, _ = lowered_problem() do model
        @variable(model, -1.0 <= x <= 1.0, start = -0.75)
        @variable(model, -1.0 <= u <= 1.0)
        @constraint(model, ∂(x) == u)
        return @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    end

    # A `final` set alone infers a reach problem; with no horizon set it is infinite.
    @test problem isa PR.OptimalControlProblem
    @test problem.time === PR.Infinity()
    @test problem.state_cost === nothing
    @test problem.transition_cost === nothing

    # `start = -0.75` becomes a *singleton* initial set, not a region.
    @test LazySets.low(problem.initial_set, 1) == -0.75
    @test LazySets.high(problem.initial_set, 1) == -0.75

    @test LazySets.low(problem.target_set, 1) == -0.5
    @test LazySets.high(problem.target_set, 1) == 0.5

    # The system is continuous-time, with X/U taken from the variable bounds.
    @test problem.system isa MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem
    @test LazySets.dim(problem.system.X) == 1
    @test LazySets.dim(problem.system.U) == 1
    @test LazySets.low(problem.system.U, 1) == -1.0
    @test LazySets.high(problem.system.U, 1) == 1.0
end

@testset "lowered problem: variable roles are inferred from dynamics alone" begin
    problem, _ = lowered_problem() do model
        @variable(model, -1.0 <= x[1:2] <= 1.0, start = 0.0)
        @variable(model, -2.0 <= u <= 2.0)
        @constraint(model, ∂(x[1]) == x[2])
        @constraint(model, ∂(x[2]) == u)
        return @constraint(model, final(x[1]) in MOI.Interval(-0.5, 0.5))
    end

    # A variable is a STATE iff it carries dynamics; everything else is an INPUT. `u` appears
    # only on a right-hand side, so it is an input without ever being declared as one.
    @test LazySets.dim(problem.system.X) == 2
    @test LazySets.dim(problem.system.U) == 1
    @test LazySets.low(problem.system.U, 1) == -2.0
end

@testset "an unused variable is an error, not a silent extra input" begin
    # FIXED-L4 (Phase 1). Before Phase 1 the unused `w` was silently classified as an input,
    # so the input space gained a dimension nobody asked for — and with an unbounded `w` the
    # input grid would have covered an infinite box.
    err = try
        lowered_problem() do model
            @variable(model, -1.0 <= x <= 1.0, start = 0.0)
            @variable(model, -1.0 <= u <= 1.0)
            @variable(model, -3.0 <= w <= 3.0)      # declared, never used
            @constraint(model, ∂(x) == u)
            return @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
        end
        nothing
    catch e
        e
    end

    @test err isa ErrorException
    # The message names the offending variable rather than an MOI index.
    @test occursin("w", err.msg)
end

@testset "lowered problem: discrete time, start sets and obstacles" begin
    problem, _ = lowered_problem() do model
        @variable(model, -1.0 <= x[1:2] <= 1.0)
        @variable(model, -0.5 <= u[1:2] <= 0.5)
        @constraint(model, Δ(x[1]) == x[1] + u[1])
        @constraint(model, Δ(x[2]) == x[2] + u[2])
        @constraint(model, start(x[1]) in MOI.Interval(-0.8, -0.7))
        @constraint(model, start(x[2]) in MOI.Interval(-0.8, -0.7))
        @constraint(model, final(x[1]) in MOI.Interval(-0.5, 0.5))
        @constraint(model, final(x[2]) in MOI.Interval(-0.5, 0.5))
        return @constraint(model, x[1:2] ∉ MOI.HyperRectangle([0.6, 0.6], [1.0, 1.0]))
    end

    # `Δ` selects the discrete-time system type.
    @test problem.system isa MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem

    # `start(...)` gives a genuine region, unlike the `start = v` keyword.
    @test LazySets.low(problem.initial_set, 1) == -0.8
    @test LazySets.high(problem.initial_set, 1) == -0.7

    # The obstacle is carved out of X, so X is no longer a plain box.
    @test !(problem.system.X isa LazySets.AbstractHyperrectangle)
    @test SVector(0.8, 0.8) ∉ problem.system.X
    @test SVector(0.0, 0.0) ∈ problem.system.X
end

@testset "input validation" begin
    # Mixing continuous and discrete dynamics is rejected.
    @test_throws ErrorException lowered_problem() do model
        @variable(model, -1.0 <= x[1:2] <= 1.0)
        @variable(model, -1.0 <= u[1:2] <= 1.0)
        @constraint(model, ∂(x[1]) == u[1])
        return @constraint(model, Δ(x[2]) == x[2] + u[2])
    end

    # A model with no dynamics at all is rejected at setup.
    @test_throws ErrorException lowered_problem() do model
        @variable(model, -1.0 <= x[1:2] <= 1.0)
        return @variable(model, -1.0 <= u[1:2] <= 1.0)
    end

    # A nonlinear constraint that is neither dynamics nor a `final`/`start` set is
    # unsupported. The wrapper `dump`s the offending function; suppress that noise.
    @test_throws ErrorException redirect_stdout(devnull) do
        return lowered_problem() do model
            @variable(model, -1.0 <= x[1:2] <= 1.0)
            @variable(model, -1.0 <= u[1:2] <= 1.0)
            @constraint(model, ∂(x[1]) == u[1])
            return @constraint(model, sin(x[2]) in MOI.Interval(0.0, 1.0))
        end
    end
end

@testset "solver attributes are forwarded to the inner optimizer" begin
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0)
    @variable(model, -1.0 <= u <= 1.0)
    @constraint(model, ∂(x) == u)

    grid = MP.GridFree(SVector(0.0), SVector(0.5))
    set_attribute(model, "state_grid", grid)
    set_attribute(model, "time_step", 0.25)
    set_attribute(model, "print_level", 0)

    # Round-trip through the wrapper onto the (hardwired) UniformGridAbstraction solver.
    @test get_attribute(model, "state_grid") === grid
    @test get_attribute(model, "time_step") == 0.25

    # DEFECT-L10 (Phase 4): the solver family is fixed in the constructor, with no way to
    # choose or auto-select one.
    @test backend(model).inner isa AB.UniformGridAbstraction.Optimizer

    # An attribute no solver recognises is rejected rather than silently dropped.
    @test_throws MOI.UnsupportedAttribute set_attribute(model, "no_such_attribute", 1)
end

@testset "options survive JuMP's caching layer" begin
    # `Model(...)` — the documented form — puts a caching layer in front of the optimizer and
    # calls `MOI.empty!` on it before copying the model in. Anything the wrapper stored in the
    # IR from an *attribute* was therefore wiped on every `optimize!`, while the same model
    # built with `direct_model` worked. Options live on the optimizer for exactly this reason.
    model = Model(Dionysos.Optimizer)
    @variable(model, -1.0 <= x <= 1.0, start = -0.75)
    @variable(model, -1.0 <= u <= 1.0)

    set_role!(x, Dionysos.STATE)
    set_attribute(model, "dynamics", (x, u) -> u)
    set_attribute(model, "time_domain", Dionysos.DISCRETE)
    @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))

    set_attribute(model, "horizon", 4.0)
    set_attribute(model, "time_step", 1.0)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(model, "print_level", 0)

    optimize!(model)

    problem = get_attribute(model, "concrete_problem")
    @test problem.time == 4                     # the horizon survived, in steps
    @test LazySets.dim(problem.system.X) == 1   # the declared role survived
    @test is_solved_and_feasible(model)         # the supplied dynamics survived
end

@testset "end-to-end: the canonical continuous reach pipeline still solves" begin
    # One genuine solve, kept tiny. `jump_frontend.jl` covers the same ground through the
    # caching backend; this pins the *direct* path used by every test above.
    x_start = -0.75

    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0, start = x_start)
    @variable(model, -1.0 <= u <= 1.0)
    @constraint(model, ∂(x) == u)
    @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))

    set_attribute(model, "time_step", 1.0)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(model, "print_level", 0)

    optimize!(model)

    concrete_controller = get_attribute(model, "concrete_controller")
    @test concrete_controller !== nothing
    @test controller_admissible(concrete_controller, SVector(x_start))

    # FIXED-L13 (Phase 2): the standard JuMP status queries answer, instead of the user
    # having to reach into solver internals to find a `success` flag.
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test result_count(model) == 1
    @test is_solved_and_feasible(model)
    @test !isempty(raw_status(model))

    # FIXED-L14 (Phase 2): one call runs the controller in closed loop, with the stopping
    # criterion taken from the specification.
    traj = Dionysos.simulate(model, SVector(x_start); nsteps = 50)
    @test !isempty(ST.states(traj))
    @test first(ST.states(traj)) == SVector(x_start)
    @test last(ST.states(traj)) ∈ get_attribute(model, "concrete_problem").target_set
end

@testset "status before optimize! and on an unsolvable model" begin
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0, start = -0.9)
    @variable(model, -1.0 <= u <= 1.0)
    @constraint(model, ∂(x) == u)
    @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))

    @test termination_status(model) == MOI.OPTIMIZE_NOT_CALLED
    @test primal_status(model) == MOI.NO_SOLUTION
    @test !is_solved_and_feasible(model)
    @test_throws ErrorException Dionysos.simulate(model, SVector(-0.9))

    # Every admissible input drives the state *away* from the target, so no controller can
    # reach it — while the abstraction itself is perfectly healthy. The status is
    # LOCALLY_INFEASIBLE, never INFEASIBLE: the abstraction is sound but not complete, so a
    # failure is a statement about this abstraction, not about the problem.
    unsolvable = direct_model(Dionysos.Optimizer())
    @variable(unsolvable, -1.0 <= y <= 1.0, start = -0.9)
    @variable(unsolvable, -1.0 <= v <= -0.5)
    @constraint(unsolvable, ∂(y) == v)
    @constraint(unsolvable, final(y) in MOI.Interval(0.5, 0.6))
    set_attribute(unsolvable, "time_step", 1.0)
    set_attribute(unsolvable, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(unsolvable, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
    set_attribute(unsolvable, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(unsolvable, "print_level", 0)
    optimize!(unsolvable)

    @test termination_status(unsolvable) == MOI.LOCALLY_INFEASIBLE
    @test primal_status(unsolvable) == MOI.NO_SOLUTION
    @test result_count(unsolvable) == 0
    @test !is_solved_and_feasible(unsolvable)
end

end # module TestMain
