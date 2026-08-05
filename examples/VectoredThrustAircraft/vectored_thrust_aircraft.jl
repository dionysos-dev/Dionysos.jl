using StaticArrays, JuMP, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "VectoredThrustAircraft",
        "vectored_thrust_aircraft.jl",
    ),
)

params = VectoredThrustAircraft.Params()

concrete_problem = VectoredThrustAircraft.optimal_control_problem(; params = params)

jacobian_bound = VectoredThrustAircraft.jacobian_bound(params)

N = 120
Δt = 0.05

x_seq = [
    SVector(
        -4.0 + 8.0 * (k - 1) / (N - 1),      # px
        -4.0 + 8.0 * (k - 1) / (N - 1),      # py
        0.35 * sin(2π * (k - 1) / (N - 1)),  # θ
        8.0 / ((N - 1) * Δt),                # vx
        8.0 / ((N - 1) * Δt),                # vy
        0.35 * cos(2π * (k - 1) / (N - 1)),  # ω
    ) for k in 1:N
]

u_seq = [
    SVector(
        1.5 * sin(2π * (k - 1) / (N - 1)),   # u₁
        0.8 * cos(2π * (k - 1) / (N - 1)),   # u₂
    ) for k in 1:(N - 1)
]

traj = ST.Trajectory(x_seq; inputs = u_seq)

system_plot! =
    VectoredThrustAircraft.system_plot!(; xlims = (-6.0, 6.0), ylims = (-6.0, 6.0))

Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1, 2),
    udims = (1, 2),
    Δt = Δt,
    fps = 5,
    # filename = "vectored_thrust_aircraft_dashboard.mp4",
    xlabel_state = "x",
    ylabel_state = "y",
    xlabel_input = "u₁",
    ylabel_input = "u₂",
    xlims_state = (-6.0, 6.0),
    ylims_state = (-6.0, 6.0),
)
