using StaticArrays, JuMP

using Symbolics
using MathOptSymbolicAD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

Ta, α, β = 18.0, 0.1, 2.0

model = Model(Dionysos.Optimizer);

@variable(model, 17.0 <= T <= 25.0, start = 18.0);
@variable(model, 0.0 <= u <= 1.0);

@mode(model, off);
@mode(model, on);

@constraint(off, ∂(T) == -α * (T - Ta))
@constraint(on, ∂(T) == -α * (T - Ta) + β * u)

@constraint(off, u == 0.0)
@constraint(on, 0.2 <= u <= 1.0)

add_transition!(model, off => on) do t
    return @constraint(t, T <= 19.0)
end

add_transition!(model, on => off) do t
    return @constraint(t, T >= 21.0)
end

comfortable = UT.box(SVector(20.5), SVector(23.0))

@constraint(off, [T] in Final(comfortable))
@constraint(on, [T] in Final(comfortable))

for m in (off, on)
    set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.2)))
    set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.2)))
    set_attribute(m, "time_step", 0.5)
    set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
    set_attribute(m, "jacobian_bound", u -> StaticArrays.SMatrix{1, 1}(-α))
end

optimize!(model);

termination_status(model)

is_solved_and_feasible(model)

abstract_system = get_attribute(model, "abstract_system");
concrete_problem = get_attribute(model, "concrete_problem");
concrete_controller = get_attribute(model, "concrete_controller");

concrete_problem.system

trajectory = Dionysos.simulate(model, (SVector(18.0), 1); nsteps = 60);

ST.states(trajectory)

ST.modes(trajectory)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
