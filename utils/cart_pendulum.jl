using StaticArrays
using Dionysos
using Plots

const ST = Dionysos.System

include("../problems/cart_pendulum.jl")

# ------------------------------------------------------------
# 1) Build problem
# ------------------------------------------------------------

params =
    CartPendulum.Params(; M = 1.0, m = 0.2, l = 1.0, J = 0.0, c = 0.1, γ = 0.05, g = 9.81)

problem = CartPendulum.optimal_control_problem(; params = params, objective = "swing_up")

# ------------------------------------------------------------
# 2) Fake trajectory for visualization test
# ------------------------------------------------------------

N = 120
Δt = 0.05

x_seq = [
    SVector(
        2.0 * sin(2π * (k - 1) / (N - 1)),          # p
        π * (k - 1) / (N - 1),                      # θ
        2.0 * cos(2π * (k - 1) / (N - 1)),          # v
        π / ((N - 1) * Δt),                         # ω
    ) for k in 1:N
]

u_seq = [SVector(
    5.0 * sin(4π * (k - 1) / (N - 1)),          # F
) for k in 1:(N - 1)]

x_traj = ST.Trajectory(x_seq)
u_traj = ST.Trajectory(u_seq)

# ------------------------------------------------------------
# 3) Dashboard animation
# ------------------------------------------------------------

system_plot! =
    CartPendulum.system_plot!(; params = params, xlims = (-3.0, 3.0), ylims = (-1.4, 1.4))

Dionysos.animate_trajectory_dashboard(
    system_plot!,
    x_traj,
    u_traj;
    xdims = (2, 4),      # θ vs ω
    udims = (1,),        # F over time
    Δt = Δt,
    fps = 10,
    # filename = "cart_pendulum_dashboard.mp4",
    title = "Cart-pendulum",
    xlabel_state = "θ [rad]",
    ylabel_state = "ω [rad/s]",
    xlabel_input = "time [s]",
    ylabel_input = "F [N]",
    xlims_state = (-0.1, π + 0.1),
    ylims_state = (-1.0, 1.0),
    ylims_input = (-5.5, 5.5),
)
