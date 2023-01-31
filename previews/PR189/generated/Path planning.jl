using StaticArrays

using Dionysos
using Dionysos.Problem
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "PathPlanning.jl"))

problem = PathPlanning.problem();

F_sys = problem.system.f;
_X_ = problem.system.X;
_U_ = problem.system.U;

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.2, 0.2, 0.2);
state_grid = DO.GridFree(x0, h);

u0 = SVector(0.0, 0.0);
h = SVector(0.3, 0.3);
input_grid = DO.GridFree(u0, h);

using JuMP
optimizer = MOI.instantiate(Abstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.optimize!(optimizer)

controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("controller"));

nstep = 100;
x0 = SVector(0.4, 0.4, 0.0);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

