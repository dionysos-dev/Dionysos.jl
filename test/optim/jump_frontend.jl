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

end # module TestMain
