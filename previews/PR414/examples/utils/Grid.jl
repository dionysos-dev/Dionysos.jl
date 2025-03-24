# # Grid
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Grid.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Grid.ipynb)
#
# In this file, we will show the different partition of the state space implemented
# - classical grid (composed of regular hyperrectangles)
# - deformed grid
# - nested classical grid
# - ellipsoidal partition based on a grid
# 
# First, let us import a few packages that are necessary to run this example.

using Dionysos
using StaticArrays
using LinearAlgebra
using Plots

# The main package [Dionysos](https://github.com/dionysos-dev/Dionysos.jl) provides most important data structures that we will need.

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain

# ### Classical grid (composed of regular hyperrectangles)
x0 = SVector(0.0, 0.0)
h = SVector(1.0 / 5, 1.0 / 5)
grid = DO.GridFree(x0, h)
rectX = UT.HyperRectangle(SVector(-2, -2), SVector(2, 2));
domainX = DO.DomainList(grid)
DO.add_set!(domainX, rectX, DO.INNER)
plot(; aspect_ratio = :equal);
plot!(domainX)

# ### Deformed grid
# We define some invertible transformation (with their inverse)
f1(x) = x
fi1(x) = x

f2(x) = x + SVector(5.0, 5.0)
fi2(x) = x - SVector(5.0, 5.0)

f3(x) = SVector(x[2] + sin(x[1]), x[1])
fi3(x) = SVector(x[2], x[1] - sin(x[2]))

f4(x) = SVector(x[1] * cos(x[2]), x[1] * sin(x[2]))
fi4(x) = SVector(sqrt(x[1] * x[1] + x[2] * x[2]), atan(x[2], x[1]))

function rotate(x, θ)
    R = @SMatrix [
        cos(θ) -sin(θ)
        sin(θ) cos(θ)
    ]
    return R * x
end
function build_f_rotation(θ; c = SVector(0.0, 0.0))
    f(x) = rotate(x - c, θ) + c
    fi(x) = rotate(x - c, -θ) + c
    return f, fi
end

function plot_deformed_grid_with_DomainList(f, fi)
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 2 * π))
    grid = DO.GridFree(SVector(0.0, 0.0), SVector(3.0, 2 * pi / 8.0))
    Dgrid = DO.DeformedGrid(grid, f, fi)
    dom = DO.DomainList(Dgrid)
    DO.add_set!(dom, X, DO.OUTER)
    plot(; aspect_ratio = :equal)
    return plot!(dom; show = true, color = :grey, efficient = false)
end

function plot_deformed_grid_with_GeneralDomain(f, fi)
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 10.0))
    obstacle = UT.HyperRectangle(SVector(10.0, 10.0), SVector(15.0, 15.0))
    hx = [2.5, 2.5]
    d = DO.RectangularObstacles(X, [obstacle])
    dom = DO.GeneralDomainList(hx; elems = d, f = f, fi = fi, fit = true)
    plot(; aspect_ratio = :equal)
    return plot!(dom; show = true, color = :grey, efficient = false)
end

rect = UT.HyperRectangle(SVector(0.0, 0.0), SVector(2.0, 2.0))
shape = UT.DeformedRectangle(rect, f2)
plot(; aspect_ratio = :equal)
plot!(rect; color = :grey, efficient = false, label = "Original")
plot!(shape; color = :red, efficient = false, label = "Deformed")

# Display some deformed Grids
plot_deformed_grid_with_DomainList(f1, fi1)
plot_deformed_grid_with_GeneralDomain(f1, fi1)

plot_deformed_grid_with_DomainList(f2, fi2)
plot_deformed_grid_with_GeneralDomain(f2, fi2)

plot_deformed_grid_with_DomainList(f3, fi3)
plot_deformed_grid_with_GeneralDomain(f3, fi3)

plot_deformed_grid_with_DomainList(f4, fi4)
plot_deformed_grid_with_GeneralDomain(f4, fi4)

f, fi = build_f_rotation(π / 3.0)
plot_deformed_grid_with_DomainList(f, fi)

# ### Nested classical grid
X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
obstacle = UT.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
hx = [3.0, 1.0] * 2.0
periodic = Int[1]
periods = [30.0, 30.0]
T0 = [0.0, 0.0]
d = DO.RectangularObstacles(X, [obstacle])
dom = DO.GeneralDomainList(
    hx;
    elems = d,
    periodic = periodic,
    periods = periods,
    T0 = T0,
    fit = true,
)
Ndomain = DO.NestedDomain(dom)
DO.cut_pos!(Ndomain, (2, 2), 1)
DO.cut_pos!(Ndomain, (2, 3), 1)
DO.cut_pos!(Ndomain, (4, 4), 2)

fig = plot(; aspect_ratio = :equal, legend = false);
plot!(Ndomain; color = :grey, efficient = false)

# ### Ellipsoidal partition based on a grid
x0 = SVector(0.0, 0.0)
n_step = 2
h = SVector(1.0 / n_step, 1.0 / n_step)
P = 0.5 * diagm((h ./ 2) .^ (-2))
rectX = UT.HyperRectangle(SVector(-2, -2), SVector(2, 2))

grid = DO.GridEllipsoidalRectangular(x0, h, P)
domain = DO.DomainList(grid)
DO.add_set!(domain, rectX, DO.OUTER)
plot(; aspect_ratio = :equal);
plot!(domain; color = :grey, opacity = 0.5, efficient = false)
