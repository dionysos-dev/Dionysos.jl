using StaticArrays, LinearAlgebra, Plots
import Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

using MathematicalSystems
import LazySets
const MS = MathematicalSystems

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "Pendulum",
        "simple_pendulum.jl",
    ),
)

params = SimplePendulum.Params(; l = 1.0, g = 9.81)

system = SimplePendulum.system(;
    params = params,
    _X_ = UT.box(SVector(-π, -5.0), SVector(π, 5.0)),
    _U_ = UT.box(SVector(-10.0), SVector(10.0)),
)

Δt = 0.1
discrete_time_system = ST.discretize_continuous_system(system, Δt)

input_set = MS.inputset(system)

pid_controller = ST.PIDControllers.PIDControllerVector(;
    Kp = SVector(15.0, 5.0),
    Ki = SVector(3.0, 0.0),
    Kd = SVector(0.0, 0.0),
    ref = ST.PIDControllers.ConstantSignal(SVector(3π / 4.0, 0.0)),
    error = ST.PIDControllers.WrapAnglePositionVelocityError(),
    dt = Δt,
    time_getter = ST.PIDControllers.ConstantTimeGetter(),
    umin = LazySets.low(input_set),
    umax = LazySets.high(input_set),
    e0 = SVector(0.0, 0.0),
)

nstep = 200
x0 = SVector(0.0, 0.0)

traj = ST.get_closed_loop_trajectory(discrete_time_system, pid_controller, x0, nstep)

fig = plot(; aspect_ratio = :equal)
plot!(system.X; color = :grey, hole_color = :black, opacity = 1.0, label = "")
plot!(traj; ms = 2.0, arrows = false)
display(fig)

# ------------------------------------------------------------
# Animation with dashboard
# ------------------------------------------------------------

system_plot! = SimplePendulum.system_plot!(; params = params)
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1, 2),      # phase plot θ vs ω
    udims = (1,),        # input over time
    Δt = Δt,
    fps = 5,
    # filename = "simple_pendulum_dashboard.mp4",
    xlabel_state = "θ [rad]",
    ylabel_state = "ω [rad/s]",
    xlabel_input = "time [s]",
    ylabel_input = "τ [Nm]",
    xlims_state = (-π, π),
    ylims_state = (-8, 8),
)
