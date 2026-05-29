using StaticArrays, LinearAlgebra, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

using MathematicalSystems
MS = MathematicalSystems

using JLD2

Ts = 0.005

Kp = 0.3
Ki = 0.1
Kd = 0.0

r = [-30.0, 0.0, 30.0, 0.0]

initial_state = [0.0, 0.0, 0.0, 0.0]
update_map = (x, y) -> x + Ts*(r-y[1:4])
output_map = (x, y) -> Tuple(Kp*(r-y[1:4]) + Ki*x - Kd*y[5:8])

pid_controller = ST.DiscreteDynamicController(
    initial_state,
    ST.PredicateDomain(_ -> true),
    update_map,
    output_map,
    false
)

JLD2.jldsave(".\\control_server\\scripts\\pid_controller.jld2"; controller = pid_controller)