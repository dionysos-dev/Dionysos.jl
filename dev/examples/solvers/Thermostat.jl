# # Example: a thermostat, as a hybrid system written in JuMP
#
# This is a **reachability problem** for a **hybrid system**: a room heated by a thermostat that
# switches between two discrete modes.
#
# The temperature $T$ follows
# ```math
# \dot{T} = -\alpha (T - T_a) \quad \text{(heater off)}, \qquad
# \dot{T} = -\alpha (T - T_a) + \beta u \quad \text{(heater on)},
# ```
# where $T_a$ is the ambient temperature and $u$ the heating power. The heater may switch on
# once the room has cooled to $19\,^\circ$C, and off once it has reached $21\,^\circ$C.
#
# What makes this different from the continuous examples is that the plant has **two modes**
# with different dynamics *and different admissible inputs* — the heater is simply off in one of
# them — plus **guarded transitions** between them. The whole model is written with ordinary
# JuMP constraints: a mode and a transition are scopes you attach constraints to, exactly as you
# attach them to a model.

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

# ## Definition of the problem

Ta, α, β = 18.0, 0.1, 2.0

model = Model(Dionysos.Optimizer);

# The state is the temperature and the input the heating power. Both bounds are model-wide; the
# modes narrow them below.
@variable(model, 17.0 <= T <= 25.0, start = 18.0);
@variable(model, 0.0 <= u <= 1.0);

# Two discrete modes.
@mode(model, off);
@mode(model, on);

# Each mode carries its own dynamics…
@constraint(off, ∂(T) == -α * (T - Ta))
@constraint(on, ∂(T) == -α * (T - Ta) + β * u)

# …and its own input set: the heater is off in one mode and throttled in the other. This is an
# ordinary bound constraint — it is *scoped* to the mode simply by being written on it.
@constraint(off, u == 0.0)
@constraint(on, 0.2 <= u <= 1.0)

# A transition is a scope too. A constraint written on it is a **guard** — which states let the
# switch happen — unless it contains `∂`/`Δ`, in which case it is the **reset map**. With no
# reset the state carries over unchanged, which is what switching a heater does.
add_transition!(model, off => on) do t
    return @constraint(t, T <= 19.0)
end

add_transition!(model, on => off) do t
    return @constraint(t, T >= 21.0)
end

# The goal: reach a comfortable band, in either mode. A specification written on a mode applies
# in that mode only, so stating it on both means "comfortable, whatever the heater is doing".
comfortable = UT.box(SVector(20.5), SVector(23.0))

@constraint(off, [T] in Final(comfortable))
@constraint(on, [T] in Final(comfortable))

# ## Definition of the abstraction
#
# Each mode is abstracted by its own sub-solver, so the discretization is set **per mode**.
# Options set on the model itself apply to every mode unless the mode overrides them.

for m in (off, on)
    set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.2)))
    set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.2)))
    set_attribute(m, "time_step", 0.5)
    set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
    set_attribute(m, "jacobian_bound", u -> StaticArrays.SMatrix{1, 1}(-α))
end

# Nothing selects a solver: the model has modes, so the hybrid abstraction is chosen for it.

optimize!(model);

# ## Results
#
# `termination_status` reports `LOCALLY_INFEASIBLE` rather than `INFEASIBLE` when no controller
# is found — the abstraction is sound but not complete, so a failure is a statement about *this*
# abstraction, and a finer grid may still succeed.

termination_status(model)

#-

is_solved_and_feasible(model)

#-

abstract_system = get_attribute(model, "abstract_system");
concrete_problem = get_attribute(model, "concrete_problem");
concrete_controller = get_attribute(model, "concrete_controller");

# The lowered system is a genuine `HybridSystems.HybridSystem` — the same object the solver
# consumes when it is built by hand.
concrete_problem.system

# The augmented state of a hybrid model is `(x, mode)`, so the simulation starts from a
# temperature *and* a mode: the room at 18 °C with the heater off.

trajectory = Dionysos.simulate(model, (SVector(18.0), 1); nsteps = 60);

# The temperature over the run, and the mode it was in at each step:

ST.states(trajectory)

#-

ST.modes(trajectory)
