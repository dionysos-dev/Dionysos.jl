using StaticArrays, Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems/single_pendulum.jl")) 

####################################################################################

g = 9.81
l = 1.0

function dynamic()
    return (x::SVector{2, Float64}, u::SVector{1, Float64}) -> begin
        return SVector{2}(
            x[2], 
            u[1] - g / l * sin(x[1]),
            )
    end
end

function jacobian_bound(u)
    return SMatrix{2, 2}(
        0.0, 1.0, 
        (g / l), 0,
        )
end

timesteps = collect(0.1:0.1:0.3) # try also 0.3, 0.5
m = length(timesteps)
# m = 5
####################################################################################

function transition_cost(q, u::Tuple{SVector{1, T}, T}) where {T}
    true_u, t = u
    return true_u[1]^2 * t
end

concrete_problem = SinglePendulum.problem(transition_cost = transition_cost)
concrete_system = concrete_problem.system

x0 = SVector(0.0, 0.0)
hx = 0.15
state_grid = DO.GridFree(x0, SVector(hx, hx))

u0 = SVector(0.0)
h = SVector(0.3 * m)
input_grid = DO.GridFree(u0, h)

using JuMP
model = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(model, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(model, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(model, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(model, MOI.RawOptimizerAttribute("time_step"), timesteps)
MOI.set(
    model,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    jacobian_bound,
)
# MOI.set(
#     optimizer,
#     MOI.RawOptimizerAttribute("approx_mode"),
#     AB.UniformGridAbstraction.GROWTH, # USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION
# )
# MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
# MOI.set(
#     optimizer,
#     MOI.RawOptimizerAttribute("automaton_constructor"),
#     (n, m) -> SY.NewIndexedAutomatonList(n, m),
# )
# MOI.set(
#     optimizer,
#     MOI.RawOptimizerAttribute("controller_constructor"),
#     () -> ST.SymbolicControllerDict(),
# )

MOI.optimize!(model)
####################################################################################

# Get the results
abstract_system = MOI.get(model, MOI.RawOptimizerAttribute("abstract_system"));
abstract_problem = MOI.get(model, MOI.RawOptimizerAttribute("abstract_problem"));
abstract_controller = MOI.get(model, MOI.RawOptimizerAttribute("abstract_controller"));
concrete_controller = MOI.get(model, MOI.RawOptimizerAttribute("concrete_controller"))
concrete_problem = MOI.get(model, MOI.RawOptimizerAttribute("concrete_problem"));
concrete_system = concrete_problem.system;
abstract_value_function = MOI.get(model, MOI.RawOptimizerAttribute("abstract_value_function"));
concrete_value_function = MOI.get(model, MOI.RawOptimizerAttribute("concrete_value_function"));

tstep = MOI.get(model, MOI.RawOptimizerAttribute("time_step"));
println("Time step used for simulation: $tstep")

println(typeof(concrete_system))

nstep = 100
function reached(x)
    if x ∈ concrete_problem.target_set
        return true
    else
        return false
    end
end

x0 = SVector(Dionysos.Utils.sample(concrete_problem.initial_set)...)

println("Value at initial state: ", concrete_value_function(x0))

control_trajectory = Dionysos.System.get_closed_loop_trajectory(
    concrete_system, #continuous time
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

using Plots

# Here we display the coordinate projection on the two first components of the state space along the trajectory.
fig = plot(; aspect_ratio = :equal);

# We display the concrete domain
plot!(concrete_system.X; color = :grey, label = "");

# We display the specifications
plot!(concrete_problem.initial_set; color = :green, opacity = 1.0, label = "Initial set");
plot!(concrete_problem.target_set; color = :red, opacity = 1.0, label = "Target set");
# plot!(abstract_system.Xdom; color = :blue, opacity = 0.1, label = "Abstract domain");
plot!(control_trajectory; markersize = 1, arrows = false)