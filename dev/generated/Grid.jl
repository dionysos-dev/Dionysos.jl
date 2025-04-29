using Dionysos
using StaticArrays, LinearAlgebra, Plots

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain

x0 = SVector(0.0, 0.0)
h = SVector(1.0 / 5, 1.0 / 5)
grid = DO.GridFree(x0, h)
rectX = UT.HyperRectangle(SVector(-2, -2), SVector(2, 2));
domainX = DO.DomainList(grid)
DO.add_set!(domainX, rectX, DO.INNER)
plot(; aspect_ratio = :equal);
plot!(domainX; efficient = false, color = :grey, label = "Grid")

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

rect = UT.HyperRectangle(SVector(0.0, 0.0), SVector(2.0, 2.0))
shape = UT.DeformedRectangle(rect, f2)
plot(; aspect_ratio = :equal)
plot!(rect; color = :grey, efficient = false, label = "Original")
plot!(shape; color = :red, efficient = false, label = "Deformed")

plot_deformed_grid_with_DomainList(f1, fi1)

plot_deformed_grid_with_DomainList(f2, fi2)

plot_deformed_grid_with_DomainList(f3, fi3)

plot_deformed_grid_with_DomainList(f4, fi4)

f, fi = build_f_rotation(π / 3.0)
plot_deformed_grid_with_DomainList(f, fi)

X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
obstacle = UT.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))

hx = SVector(6.0, 2.0)
periodic_dims = SVector(1)
periods = SVector(30.0)
start = SVector(0.0)

free_space = UT.LazySetMinus(X, UT.LazyUnionSetArray([obstacle]))

domain = DO.PeriodicDomainList(periodic_dims, periods, start, hx)

DO.add_set!(domain, free_space, DO.INNER)

Ndomain = DO.NestedDomain(domain)

div = 3
DO.cut_pos!(Ndomain, (2, 2), 1; div = div)
DO.cut_pos!(Ndomain, (2, 3), 1; div = div)
DO.cut_pos!(Ndomain, (4, 4), 2; div = div)

fig = plot(; aspect_ratio = :equal, legend = false)
plot!(free_space; color = :blue, label = "Base domain", efficient = false)
plot!(Ndomain; color = :grey, efficient = false)

x0 = SVector(0.0, 0.0)
n_step = 2
h = SVector(1.0 / n_step, 1.0 / n_step)
P = 0.5 * diagm((h ./ 2) .^ (-2))

rectX = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))

ellip_grid = DO.GridEllipsoidalRectangular(x0, h, P)
ellip_domain = DO.DomainList(ellip_grid)
DO.add_set!(ellip_domain, rectX, DO.OUTER)

fig2 = plot(; aspect_ratio = :equal)
plot!(ellip_domain; color = :grey, opacity = 0.5, efficient = false)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
