using StaticArrays

using Dionysos
using Dionysos.Problem
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "DCDC.jl"))

problem = DCDC.problem()

x0 = SVector(0.0, 0.0)
hx = SVector(2.0/4.0e3, 2.0/4.0e3)
state_grid = DO.GridFree(x0, hx)
u0 = SVector(1)
hu = SVector(1)
input_grid = DO.GridFree(u0, hu)

using JuMP
optimizer = MOI.instantiate(Abstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.optimize!(optimizer)

controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("controller"));

#nstep = 300;
#x0 = SVector(1.2, 5.6);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

