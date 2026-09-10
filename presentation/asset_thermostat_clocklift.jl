# Lifts 2/4: the SAME 1-D thermostat plant as slide 1, abstracted with the
# reusable ClockLift: the solver's `clock` attribute wraps the spatial
# abstraction into an (T, t) product where the clock advances exactly.
# The task is a timed one: reach the band 21-23 °C with the clock in [3, 4] s.
#
# The concrete closed loop over a clock-lifted continuous plant is not wired in
# the engine yet, so the rollout below queries the controller at (T, t) and
# integrates the affine plant exactly.
#
# Output: presentation/gallery/thermostat_clocklift.gif
# Run: julia --project=test presentation/asset_thermostat_clocklift.jl

ENV["GKSwstype"] = "100"

using StaticArrays, Plots
import LazySets
import MathOptInterface as MOI
using Dionysos
const DI = Dionysos
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction

include(joinpath(dirname(@__DIR__), "problems", "Thermostat", "thermostat_system.jl"))

const Δt = 0.1
const t_max = 5.0
const min_time, max_time = 3.0, 4.0

system = ThermostatSystem.system()          # 1-D plant, T ∈ [18, 24]
initial = LazySets.Hyperrectangle(; low = [18.0, 0.0], high = [18.5, 0.0])
target = LazySets.Hyperrectangle(; low = [21.0, min_time], high = [23.0, max_time])
problem = PR.OptimalControlProblem(system, initial, target, nothing, nothing, PR.Infinity())

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("state_grid"),
    MP.GridFree(SVector(0.0), SVector(0.05)),
)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("input_grid"),
    MP.GridFree(SVector(1.0), SVector(1.0)),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    ThermostatSystem.jacobian_bound(),
)
# The reusable time lift: the clock is a factor of the abstraction, not a state
# dimension, so it advances exactly one slice per step.
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("clock"),
    SY.ClockAbstraction(LazySets.Hyperrectangle(; low = [0.0], high = [t_max]), Δt),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
println("clock-lifted model: ", typeof(abstract_system))
println("states: ", SY.get_n_state(abstract_system))
println("base states: ", SY.get_n_state(abstract_system.base))
println("transitions: ", SY.get_n_transitions(abstract_system))
println("base transitions: ", SY.get_n_transitions(abstract_system.base))
println("success: ", MOI.get(optimizer, MOI.RawOptimizerAttribute("success")))

# ---- closed loop over (T, t), integrating the affine plant exactly ----------
p = ThermostatSystem.Params()
step_T(T, on) = begin
    b = on ? p.alpha * p.Ta + p.beta : p.alpha * p.Ta
    e = exp(-p.alpha * Δt)
    return e * T + (1 - e) / p.alpha * b
end

function rollout()
    T, t = 18.2, 0.0
    entered = false
    Ts, us, ts = [SVector(T)], SVector{1, Float64}[], [t]
    for _ in 1:Int(round(t_max / Δt))
        u = ST.output_control(controller, nothing, SVector(T, t))
        if u === nothing
            println("controller undefined at t = ", round(t; digits = 2), " s")
            break
        end
        on = Float64(u[1]) == 1.0
        T = step_T(T, on)
        t += Δt
        push!(us, SVector(Float64(u[1])))
        push!(Ts, SVector(T))
        push!(ts, t)
        if 21.0 <= T <= 23.0 && min_time <= t <= max_time && !entered
            entered = true
            println(
                "entered the band at t = ",
                round(t; digits = 2),
                " s, T = ",
                round(T; digits = 2),
                " °C",
            )
        end
    end
    return Ts, us, ts
end

Ts, us, ts = rollout()

trajectory = ST.Trajectory(Ts; inputs = us, times = ts)

system_plot! = ThermostatSystem.system_plot!(; problem = problem)
anim = DI.animate_trajectory_dashboard(
    system_plot!,
    trajectory;
    xdims = (1,),
    udims = (1,),
    Δt = Δt,
    frame_step = 2,
    title = "Clock-lifted thermostat (T, t)",
    ylabel_state = "T [°C]",
    ylabel_input = "u (1 = ON, 2 = OFF)",
)
out = joinpath(@__DIR__, "gallery", "thermostat_clocklift.gif")
gif(anim, out; fps = 6)
println("saved ", out)
