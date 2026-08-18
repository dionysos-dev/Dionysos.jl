using StaticArrays
using JuMP
using Symbolics, MathOptSymbolicAD # compile the ∂ dynamics of the front-end
using Dionysos
import LazySets
import Plots

const DI = Dionysos
const MP = DI.Mapping
const ST = DI.System

const L1 = 0.20125 # hip → knee
const L2 = 0.172   # knee → foot

function joint_positions(θ::SVector{4})
    knee_s = L2 * SVector(-sin(θ[1] + θ[2]), cos(θ[1] + θ[2]))     # stance knee
    hip = knee_s + L1 * SVector(-sin(θ[1]), cos(θ[1]))
    knee_w = hip + L1 * SVector(sin(θ[3]), -cos(θ[3]))             # swing knee
    foot = knee_w + L2 * SVector(sin(θ[3] + θ[4]), -cos(θ[3] + θ[4]))
    return (; knee_s, hip, knee_w, foot)
end

swing_foot(θ::SVector{4}) = joint_positions(θ).foot

dx, τ = 0.1, 0.1
du = dx / τ                       # one speed notch; u ∈ {-du, 0, du} per joint
state_grid = MP.GridFree(SVector(0.0, 0.0, 0.0, 0.0), SVector(dx, dx, dx, dx))
input_grid = MP.GridFree(SVector(0.0, 0.0, 0.0, 0.0), SVector(du, du, du, du))
MP.is_lattice_exact(state_grid, input_grid, τ)

obstacle =
    LazySets.VPolygon([SVector(-0.39, 0.0), SVector(-0.35, 0.02), SVector(-0.31, 0.0)])

x_bar = 0.6
X_box = LazySets.Hyperrectangle(;
    low = SVector(-x_bar, -x_bar, -x_bar, -x_bar),
    high = SVector(x_bar, x_bar, x_bar, x_bar),
)

grad = SVector(L1 + L2, L2, L1 + L2, L2)          # per-joint Lipschitz bound
dev = sum(grad .* MP.get_h(state_grid) ./ 2)

# `MP.cells_where` collects grid cells by predicate into a `MP.CellUnion` — a
# cell-aligned set whose discretization is exact (a hole removes exactly its
# cells, a target is recovered exactly). The obstacle test could also use
# `MP.image_blocked_cells(swing_foot, grad, ...)` directly; it is spelled out
# here to keep every ingredient visible.
removed = MP.cells_where(state_grid, X_box) do pos
    foot = swing_foot(SVector{4}(MP.get_coord_by_pos(state_grid, pos)))
    below_ground = foot[2] < -dev
    hits_obstacle =
        !Dionysos.Utils.is_disjoint(
            LazySets.Hyperrectangle(foot, SVector(dev, dev)),
            obstacle,
        )
    return below_ground || hits_obstacle
end
length(removed)

foothold = SVector(0.2, 0.0)

target_cells = Set{NTuple{4, Int}}()
for θ1 in range(-x_bar, x_bar; step = dx / 3), θ2 in range(-x_bar, x_bar; step = dx / 3)
    hip = joint_positions(SVector(θ1, θ2, 0.0, 0.0)).hip
    d = foothold - hip
    c4 = (d[1]^2 + d[2]^2 - L1^2 - L2^2) / (2 * L1 * L2)
    -1.0 <= c4 <= 1.0 || continue
    θ4 = -acos(c4)
    θ3 = atan(d[1], -d[2]) - atan(L2 * sin(θ4), L1 + L2 * cos(θ4))
    θ = SVector(θ1, θ2, θ3, θ4)
    all(abs.(θ) .<= x_bar) || continue
    push!(target_cells, MP.get_pos_by_coord(state_grid, θ))
end
target_set = MP.CellUnion(state_grid, target_cells)
length(target_cells)

x0 = SVector(0.2, 0.0, -0.2, 0.0)

model = Model(Dionysos.Optimizer)
@variable(model, -x_bar <= x[i = 1:4] <= x_bar, start = x0[i])
@variable(model, -du <= u[1:4] <= du)
@constraint(model, [i = 1:4], ∂(x[i]) == u[i])
@constraint(model, x ∉ removed)
@constraint(model, x in Final(target_set))

set_attribute(model, "time_step", τ)
set_attribute(model, "state_grid", state_grid)
set_attribute(model, "input_grid", input_grid)
set_attribute(
    model,
    "approx_mode",
    DI.Optim.Abstraction.UniformGridAbstraction.CENTER_SIMULATION,
)

optimize!(model)
get_attribute(model, "success")

controller_free = get_attribute(model, "concrete_controller")
discrete_time_system = get_attribute(model, "discrete_time_system")

reached(xk) = xk ∈ target_set
traj_free = ST.get_closed_loop_trajectory(
    discrete_time_system,
    controller_free,
    x0,
    100;
    stopping = reached,
)
xs_free = collect(ST.states(traj_free))
(length(xs_free) - 1, reached(xs_free[end]))

rest = SVector(0.0, 0.0, 0.0, 0.0)
slew = DI.Optim.DiscreteSystems.BoundedInputVariation(
    (u1, u2) -> maximum(abs.(u1 - u2)),
    du;
    target_input = rest,
    initial_input = rest,
)
set_attribute(model, "bounded_input_variation", slew)
optimize!(model)
get_attribute(model, "success")

controller = get_attribute(model, "concrete_controller")
traj = ST.get_closed_loop_trajectory(
    discrete_time_system,
    controller,
    x0,
    100;
    stopping = reached,
)
xs = collect(ST.states(traj))
us = collect(ST.inputs(traj))
(length(xs) - 1, reached(xs[end]))

maximum(maximum(abs.(us[k + 1] - us[k])) for k in 1:(length(us) - 1))

function draw_robot!(fig, θ::SVector{4})
    p = joint_positions(θ)
    Plots.plot!(
        fig,
        [0.0, p.knee_s[1], p.hip[1]],
        [0.0, p.knee_s[2], p.hip[2]];
        lw = 4,
        marker = :circle,
        color = :steelblue,
        label = "",
    )
    Plots.plot!(
        fig,
        [p.hip[1], p.knee_w[1], p.foot[1]],
        [p.hip[2], p.knee_w[2], p.foot[2]];
        lw = 4,
        marker = :circle,
        color = :indianred,
        label = "",
    )
    return fig
end

fig = Plots.plot(; aspect_ratio = :equal, legend = false, xlims = (-0.6, 0.6))
Plots.plot!(fig, [-0.6, 0.6], [0.0, 0.0]; color = :black, lw = 1, label = "")
Plots.plot!(fig, obstacle; color = :black, alpha = 0.8, label = "")
Plots.scatter!(fig, [foothold[1]], [foothold[2]]; marker = :xcross, color = :green)
for (k, xk) in enumerate(xs)
    (k % 4 == 1 || k == length(xs)) && draw_robot!(fig, xk)
end
fig

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
