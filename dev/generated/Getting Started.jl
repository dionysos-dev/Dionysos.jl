using Dionysos
using StaticArrays, MathematicalSystems
using LinearAlgebra
using Plots

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic

rectX = UT.HyperRectangle(SVector(-2, -2), SVector(2, 2));
rectU = UT.HyperRectangle(SVector(-5), SVector(5));

x0 = SVector(0.0, 0.0)
hx = SVector(1.0/5, 1.0/5)
Xgrid = MP.GridFree(x0, hx);

u0 = SVector(0.0)
hu = SVector(1.0/5)
Ugrid = MP.GridFree(u0, hu);

Xmap = MP.ExplicitGridMapping(Xgrid)
Umap = MP.ExplicitGridMapping(Ugrid)

MP.add_set!(Xmap, rectX, MP.INNER)
MP.add_set!(Umap, rectU, MP.INNER)

Xset = MP.MappingSet{2}()  # default "all states"
Uset = MP.MappingSet{1}()

tstep = 0.1

A = SMatrix{2, 2}(0.0, 1.0, -3.0, 1.0)
B = SMatrix{2, 1}(0.0, 1.0)

F_sys = (x, u) -> A*x + B*u

jacobian_bound = u -> abs.(A)

concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
    F_sys,
    2,
    1,
    nothing,
    nothing,
)

continuous_approx =
    ST.ContinuousTimeGrowthBound_from_jacobian_bound(concrete_system, jacobian_bound)
discrete_approx = ST.discretize(continuous_approx, tstep)

concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
    F_sys,
    2,
    1,
    nothing,
    nothing,
)
continuous_approx =
    ST.ContinuousTimeGrowthBound_from_jacobian_bound(concrete_system, jacobian_bound)
discrete_approx = ST.discretize(continuous_approx, tstep)

abstract_system = SY.SymbolicModelList(
    Xmap,
    Umap;
    Xset = Xset,
    Rset = Xset,  # allowed targets
    Uset = Uset,
)
SY.compute_abstract_system_from_concrete_system!(abstract_system, discrete_approx)

xpos = MP.get_pos_by_coord(Xgrid, SVector(1.1, 1.3))
x_center = MP.get_coord_by_pos(Xgrid, xpos)

q = MP.get_state_by_pos(Xmap, xpos)
abstract_input = 1
u = SY.get_concrete_input(abstract_system, abstract_input)

post = Int[]
ST.compute_post!(post, SY.get_automaton(abstract_system), q, abstract_input)

post_set = SY.get_state_set_from_states(abstract_system, post)

fig = plot(; aspect_ratio = :equal)
dims = [1, 2]
xlims!(-2, 2);
ylims!(-2, 2)

plot!((Xset, Xmap); fc = "grey", dims = dims, label = "X", efficient = false)
plot!(
    (MP.stateset_from_states(Xmap, [q]), Xmap);
    fc = "blue",
    dims = dims,
    label = "x cell",
)
plot!((post_set, Xmap); fc = "green", dims = dims, label = "Post", efficient = false)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
