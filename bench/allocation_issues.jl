using Dionysos
using MathematicalSystems
using StaticArrays
using LinearAlgebra

const DO = Dionysos.Domain
const UT = Dionysos.Utils
const ST = Dionysos.System
const MS = MathematicalSystems
const SY = Dionysos.Symbolic

# ==========================================
# =========[CREATION OF THE SYSTEM]=========
# ==========================================
# State domain [0,1]^2
lb = SVector(0.0, 0.0)
ub = SVector(1.0, 1.0)
h = (ub - lb) ./ (10 - 1)
Xgrid = DO.GridFree(lb, h)
Xfull = DO.DomainList(Xgrid)
DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)
# Input domain [-1,1]^2 
lb_u = SVector(-1.0, -1.0)
ub_u = SVector(1.0, 1.0)
h_u = SVector(0.5, 0.5)
Ugrid = DO.GridFree(lb_u, h_u)
Ufull = DO.DomainList(Ugrid)
DO.add_set!(Ufull, UT.HyperRectangle(lb_u, ub_u), DO.OUTER)

A = @SMatrix [
    0.0 1.0;
    1.0 0.0;
]
B = @SMatrix [
    1.0 0.0;
    0.0 1.0;
]

# version 1 ; numb of alloc : 6 151 919
# function F_dyn(x::SVector{2,Float64}, u::SVector{2,Float64})
#     @inbounds return SVector{2,Float64}(
#         A[1,1]*x[1] + A[1,2]*x[2] + B[1,1]*u[1] + B[1,2]*u[2],
#         A[2,1]*x[1] + A[2,2]*x[2] + B[2,1]*u[1] + B[2,2]*u[2]
#     )
# end 

# version 2 ; num of alloc = 5 124 299
# F_dyn(x, u) = A * x + B * u  

# version 3 ; num of alloc = 6 341 417
# F_dyn(x, u) = SVector(
#     A[1,1] * x[1] + A[1,2] * x[2] + B[1,1] * u[1] + B[1,2] * u[2],
#     A[2,1] * x[1] + A[2,2] * x[2] + B[2,1] * u[1] + B[2,2] * u[2]
# ) 

# version 4 ; numb of alloc : 4 644 861 
F_dyn(x::SVector{2,Float64}, u::SVector{2,Float64}) = @inbounds SVector(x[2] + u[1], x[1] + u[2]) 


continuous_system =
    MS.ConstrainedBlackBoxControlContinuousSystem(F_dyn, 2, 2, nothing, nothing)
cont_center = ST.ContinuousTimeCenteredSimulation(continuous_system)
centered = ST.discretize(cont_center, 0.5)
symmodel_builder = () -> SY.NewSymbolicModelListList(Xfull, Ufull)
abstract_sys = symmodel_builder()

# ==========================================
# ========[solve allocation issues]=========
# ==========================================
# println("\n=== Print number of alloc in symmodel_builder() ===")
# t = @timed abstract_sys = symmodel_builder()
# println("numb of alloc : ", t.gcstats.poolalloc)

# println("\n=== Profiling compute_abstract_system_from_concrete_system ===")
# using Profile
# Profile.clear()
# @profile SY.compute_abstract_system_from_concrete_system!(abstract_sys, centered)
# Profile.print()


# println("\n=== Print number of alloc in compute_abstract_system_from_concrete_system ===")
# t = @timed _ = SY.compute_abstract_system_from_concrete_system!(abstract_sys, centered)
# println("numb of alloc : ", t.gcstats.poolalloc)


println("\n\navant warm-up:")
GC.gc() # Force garbage collection to get a clean state
t = @timed _ = SY.compute_abstract_system_from_concrete_system!(abstract_sys, centered)
println("\n numb of alloc : ", t.gcstats.poolalloc)
GC.gc() # Force garbage collection to get a clean state
println("\n\naprès warm-up:")
t = @timed _ = SY.compute_abstract_system_from_concrete_system!(abstract_sys, centered)
println("numb of alloc : ", t.gcstats.poolalloc)


