using StaticArrays, Plots
using Dionysos
import LazySets
const DI = Dionysos
const UT = DI.Utils
const MP = DI.Mapping

function visualize_explicit_mapping(m::MP.PeriodicGridMapping)
    grid = MP.get_grid(m)

    p = SVector(-1.5, 5.5)

    pos_raw = MP.get_pos_by_coord(grid, p)
    coord = MP.wrap_coord(m, p)

    # Add a rectangle that spans all x (periodic) at fixed y = 0
    rect = LazySets.Hyperrectangle(; low = SVector(4.5, 2.5), high = SVector(4.8, 4.5))
    MP.cover!(m, rect, MP.OUTER)

    fig = plot(; aspect_ratio = :equal)
    plot!(m; label = "", color = :grey, efficient = false)
    plot!(rect; label = "", color = :blue, opacity = 0.5, efficient = false)
    plot!(MP.get_elem_by_coord(m, p); color = :green, label = "", efficient = false)
    scatter!(
        [p[1]],
        [p[2]];
        label = "original coord p",
        color = :red,
        marker = :cross,
        markersize = 8,
    )
    scatter!(
        [coord[1]],
        [coord[2]];
        label = "wrapped coord",
        color = :green,
        marker = :star5,
        markersize = 10,
    )
    return display(fig)
end

function visualize_implicit_mapping(m::MP.PeriodicGridMapping)
    grid = MP.get_grid(m)

    p = SVector(-1.5, 5.5)

    pos_raw = MP.get_pos_by_coord(grid, p)
    coord = MP.wrap_coord(m, p)

    # # Add a rectangle that spans all x (periodic) at fixed y = 0
    # rect = LazySets.Hyperrectangle(;
    #     low = SVector(4.5, 2.5),
    #     high = SVector(4.8, 4.5),
    # )
    # MP.add_set!(m, rect, MP.OUTER)

    fig = plot(; aspect_ratio = :equal)
    plot!(m; label = "", color = :grey, efficient = false)
    plot!(MP.get_elem_by_coord(m, p); color = :green, label = "", efficient = false)
    scatter!(
        [p[1]],
        [p[2]];
        label = "original coord p",
        color = :red,
        marker = :cross,
        markersize = 8,
    )
    scatter!(
        [coord[1]],
        [coord[2]];
        label = "wrapped coord",
        color = :green,
        marker = :star5,
        markersize = 10,
    )
    return display(fig)
end

# -----------------------------------------------------
# ExplicitGridMapping
# -----------------------------------------------------

# ----------------------------
# Domain setup: 2D with dim 1 periodic
# ----------------------------
# periodic_dims = SVector(1)
# periods = SVector(4.0)
# start = SVector(0.0)

# # Consruct the grid yourself
# origin = SVector(0.5, 0.0)
# h = SVector(1.0, 2.0)
# grid = MP.GridFree(origin, h)

# # Consruct the grid with some periodic helper
# grid = MP.get_grid_in_periods(periodic_dims, periods, start, h)

# # Construct the mapping
# mapping = MP.ExplicitGridMapping(grid)
# m = MP.PeriodicGridMapping(periodic_dims, periods, start, mapping)
# visualize_explicit_mapping(m)

# ----------------------------
# Domain setup: 2D with no dim periodic
# ----------------------------
# m = MP.PeriodicGridMapping(mapping)
# visualize_explicit_mapping(m)

# ----------------------------
# Domain setup: 2D with both dimensions periodic
# ----------------------------
# periodic_dims = SVector(1, 2)
# periods = SVector(4.0, 3.0)
# h = SVector(1.0, 1.5)
# grid = MP.get_grid_in_periods(periodic_dims, periods, h)
# mapping = MP.ExplicitGridMapping(grid)
# m = MP.PeriodicGridMapping(periodic_dims, periods, mapping)
# visualize_explicit_mapping(m)

# -----------------------------------------------------
# ImplicitGridMapping
# -----------------------------------------------------

periodic_dims = SVector(1)
periods = SVector(4.0)
start = SVector(0.0)

# Consruct the grid yourself
origin = SVector(0.5, 0.0)
h = SVector(1.0, 2.0)
grid = MP.GridFree(origin, h)

# Consruct the grid with some periodic helper
grid = MP.get_grid_in_periods(periodic_dims, periods, start, h)

# Construct the mapping
rect = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(4.0, 8.0))
mapping = MP.ImplicitGridMapping(grid, rect; incl_mode = MP.INNER)
m = MP.PeriodicGridMapping(periodic_dims, periods, start, mapping)
visualize_implicit_mapping(m)
