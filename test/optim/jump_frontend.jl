module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

# The canonical user entry point: a JuMP model on Dionysos.Optimizer, with dynamics
# expressed via `∂` and a target via `final`. `using Dionysos` brings the exported
# operators (∂, final) into scope for the @constraint macros. `Dionysos.Optimizer`
# itself is provided by the MathOptSymbolicAD extension, so those must be loaded.
using Dionysos
using JuMP
using Symbolics
using MathOptSymbolicAD
import MathOptInterface as MOI

@testset "JuMP frontend: reach on a 2D single integrator" begin
    # ẋ = u on [-1, 1]^2, inputs in [-1, 1]^2. Reach a box around the origin.
    x_start = [-0.75, -0.75]

    model = Model(Dionysos.Optimizer)
    @variable(model, -1.0 <= x[i = 1:2] <= 1.0, start = x_start[i])
    @variable(model, -1.0 <= u[1:2] <= 1.0)

    @constraint(model, ∂(x[1]) == u[1])
    @constraint(model, ∂(x[2]) == u[2])

    @constraint(model, final(x[1]) in MOI.Interval(-0.5, 0.5))
    @constraint(model, final(x[2]) in MOI.Interval(-0.5, 0.5))

    set_attribute(model, "time_step", 1.0)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.25, 0.25)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))

    optimize!(model)

    concrete_controller = get_attribute(model, "concrete_controller")
    abstract_controller = get_attribute(model, "abstract_controller")

    @test concrete_controller !== nothing
    @test abstract_controller !== nothing
    # The synthesized controller must act at the declared start state.
    @test controller_admissible(concrete_controller, SVector(x_start...))
end

@testset "JuMP frontend: discrete-time reach with start() and an obstacle" begin
    # Discrete-time counterpart of the continuous test above: dynamics via `Δ`
    # (xₖ₊₁ = xₖ + u), the initial set via `start(...)`, and a corner obstacle via `∉`.
    # Together these exercise the discrete branch of the MOI wrapper, the `:start`
    # constraint head, and the obstacle (`OuterSet`) path.
    x_start = [-0.75, -0.75]

    model = Model(Dionysos.Optimizer)
    @variable(model, -1.0 <= x[1:2] <= 1.0)
    @variable(model, -0.5 <= u[1:2] <= 0.5)

    @constraint(model, Δ(x[1]) == x[1] + u[1])
    @constraint(model, Δ(x[2]) == x[2] + u[2])

    @constraint(model, start(x[1]) in MOI.Interval(x_start[1] - 0.1, x_start[1] + 0.1))
    @constraint(model, start(x[2]) in MOI.Interval(x_start[2] - 0.1, x_start[2] + 0.1))

    @constraint(model, final(x[1]) in MOI.Interval(-0.5, 0.5))
    @constraint(model, final(x[2]) in MOI.Interval(-0.5, 0.5))

    # Obstacle in the top-right corner, clear of both the start and the target.
    @constraint(model, x[1:2] ∉ MOI.HyperRectangle([0.6, 0.6], [1.0, 1.0]))

    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.25, 0.25)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))

    optimize!(model)

    concrete_problem = get_attribute(model, "concrete_problem")
    concrete_controller = get_attribute(model, "concrete_controller")

    @test concrete_controller !== nothing
    # The obstacle was carved out of the state domain (set-minus), so it is not covered.
    @test SVector(0.8, 0.8) ∉ concrete_problem.system.X
    @test controller_admissible(concrete_controller, SVector(x_start...))
end

@testset "JuMP frontend: input validation" begin
    # Mixing continuous `∂` and discrete `Δ` dynamics on the same model is rejected.
    model = Model(Dionysos.Optimizer)
    @variable(model, -1.0 <= x[1:2] <= 1.0)
    @variable(model, -1.0 <= u[1:2] <= 1.0)
    @constraint(model, ∂(x[1]) == u[1])
    @constraint(model, Δ(x[2]) == x[2] + u[2])
    @test_throws ErrorException optimize!(model)

    # A model with no dynamics at all is rejected at setup.
    model = Model(Dionysos.Optimizer)
    @variable(model, -1.0 <= x[1:2] <= 1.0)
    @variable(model, -1.0 <= u[1:2] <= 1.0)
    @test_throws ErrorException optimize!(model)

    # A nonlinear constraint that is neither dynamics nor a `final`/`start` set is
    # unsupported. The wrapper `dump`s the offending function; suppress that noise.
    model = Model(Dionysos.Optimizer)
    @variable(model, -1.0 <= x[1:2] <= 1.0)
    @variable(model, -1.0 <= u[1:2] <= 1.0)
    @constraint(model, ∂(x[1]) == u[1])
    @constraint(model, sin(x[2]) in MOI.Interval(0.0, 1.0))
    @test_throws ErrorException redirect_stdout(devnull) do
        return optimize!(model)
    end
end

end # module TestMain
