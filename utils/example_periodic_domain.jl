using StaticArrays, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain

function visualize(domain)
    p = SVector(-1.5, 5.5)
    pos = DO.get_pos_by_coord(domain, p)
    coord = DO.wrap_coord(domain, p)

    # Add a rectangle that spans all x (periodic) at fixed y = 0
    rect = UT.HyperRectangle(
        SVector(4.5, 2.5),  # lower bound
        SVector(4.8, 3.5),    # upper bound → covers x = [0,4], y = [−0.5, 0.5]
    )
    DO.add_set!(domain, rect, DO.OUTER)

    fig = plot(; aspect_ratio = :equal);
    plot!(domain; label = "", color = :grey, efficient = false);
    plot!(rect; label = "", color = :blue, opacity = 0.5, efficient = false);
    plot!(DO.get_elem_by_coord(domain, p); color = :green, label = "", efficient = false)
    scatter!(
        [p[1]],
        [p[2]];
        label = "original coord p",
        color = :red,
        marker = :cross,
        markersize = 8,
    );
    scatter!(
        [coord[1]],
        [coord[2]];
        label = "wrapped coord",
        color = :green,
        marker = :star5,
        markersize = 10,
    );
    return display(fig)
end

# ----------------------------
# Domain setup: 2D with dim 1 periodic
# ----------------------------
periodic_dims = SVector(1)
periods = SVector(4.0)
h = SVector(1.0, 1.0)
domain = DO.PeriodicDomainList(periodic_dims, periods, h)
visualize(domain)

# ----------------------------
# Domain setup: 2D with no dim periodic
# ----------------------------
domain = DO.PeriodicDomainList(h)
visualize(domain)

# ----------------------------
# Domain setup: 2D with both dimensions periodic
# ----------------------------
periodic_dims = SVector(1, 2)
periods = SVector(4.0, 3.0)
h = SVector(1.0, 1.0)
domain = DO.PeriodicDomainList(periodic_dims, periods, h)
visualize(domain)
