using StaticArrays, LinearAlgebra, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

using MathematicalSystems
MS = MathematicalSystems

using JLD2

tstep = 0.01

pid_controller = ST.PIDControllers.PIDControllerVector(;
    Kp = 0.6649,
    Ki = 0.2817,
    Kd = 0.3923,
    ref = ST.PIDControllers.ConstantSignal(2.0),
    error = (x, r, t) -> r .- x[1],
    e0 = 0.0,
    umin = -10.0,
    umax = 10.0,
    dt = tstep,
    time_getter = ST.PIDControllers.ConstantTimeGetter()
)

JLD2.jldsave(".\\control_server\\scripts\\pid_controller.jld2"; controller = pid_controller)