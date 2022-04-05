module TestMain

using Test
using StaticArrays
using ..Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "GridDomain" begin
orig = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
grid = DO.GridFree(orig, h)
domain1 = DO.DomainList(grid)

DO.add_coord!(domain1, SVector(1.2, 3.5))
@test DO.get_ncells(domain1) == 1
DO.add_coord!(domain1, SVector(-0.5002, -1.0))
@test DO.get_ncells(domain1) == 2

DO.add_set!(domain1, UT.HyperRectangle(SVector(1.0, 0.0), SVector(11.0, 10.0)), DO.OUTER)
@test DO.get_ncells(domain1) == 67

DO.remove_coord!(domain1, SVector(2.0, 2.0))
@test DO.get_ncells(domain1) == 66
DO.remove_set!(domain1, UT.HyperRectangle(SVector(5.0, 5.0), SVector(10000.0, 10000.0)), DO.INNER)
@test DO.get_ncells(domain1) == 48

pos_iter = DO.enum_pos(domain1)
@test length(pos_iter) == 48

domain2 = DO.DomainList(grid)
union!(domain2, domain1)
@test DO.get_ncells(domain2) == 48
DO.remove_set!(domain2, UT.HyperRectangle(SVector(1.0, 1.0), SVector(2.0, 2.0)), DO.OUTER)
@test DO.get_ncells(domain2) == 45
DO.add_subset!(domain2, domain1, UT.HyperRectangle(SVector(0.0, 0.0), SVector(5.0, 5.0)), DO.INNER)
@test DO.get_ncells(domain2) == 46

# @static if get(ENV, "CI", "false") == "false"
#     include("../../src/Abstraction/plotting.jl")
#     using PyPlot
#     fig = PyPlot.figure()
#     ax = fig.gca()
#     ax.set_xlim((-2.0, 14.0))
#     ax.set_ylim((-2.0, 14.0))
#     Plot.domain!(ax, 1:2, domain2, fa = 0.3)
#     Plot.set!(ax, 1:2, AB.HyperRectangle(SVector(1.0, 0.0), SVector(8.0, 10.0)))
#     Plot.set!(ax, 1:2, AB.HyperRectangle(SVector(5.0, 5.0), SVector(10000.0, 10000.0)))
#     Plot.set!(ax, 1:2, AB.HyperRectangle(SVector(1.0, 1.0), SVector(2.0, 2.0)))
#     Plot.set!(ax, 1:2, AB.HyperRectangle(SVector(0.0, 0.0), SVector(5.0, 5.0)))
# end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
