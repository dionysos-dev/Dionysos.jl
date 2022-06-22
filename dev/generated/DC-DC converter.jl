using StaticArrays

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "DCDC.jl"))

contsys=DCDC.system()

#The default values are included, but can be changed.
DCDC.solveproblem(contsys)

#nstep = 300;
#x0 = SVector(1.2, 5.6);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

