using MathematicalSystems, StaticArrays, Plots
using JuMP
using JLD2

# include Dionysos
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems/Biped_robot/Simplified_problem",
        "robot_problem.jl",
    ),
)

#######################################################
################### File Parameters ###################
#######################################################
filename_save = joinpath(@__DIR__, "Abstraction_solver.jld2")
do_empty_optim = true
verify_save = false

#######################################################
################### Optim Parameters ##################
#######################################################
concrete_problem = RobotProblem.problem(; tstep = 1e-1)
concrete_system = concrete_problem.system

### Set the optimizer
n_state = MathematicalSystems.statedim(concrete_system)
n_input = MathematicalSystems.inputdim(concrete_system)
state_space = MathematicalSystems.stateset(concrete_system)
input_space = MathematicalSystems.inputset(concrete_system)
println("n_state: ", n_state)
println("n_input: ", n_input)

x0 = SVector{n_state, Float64}(zeros(n_state))
hx = SVector{n_state, Float64}([fill(3*π/180, 3)..., fill(0.15, 3)...])
state_grid = DO.GridFree(x0, hx)

u0 = SVector{n_input, Float64}(zeros(n_input))
hu = SVector{n_input, Float64}(fill(1.0, n_input))
input_grid = DO.GridFree(u0, hu)

using JuMP
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.Silent(), true)  
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)


### Optimize
if(do_empty_optim)
    MOI.optimize!(optimizer)

    # Save the abstraction solver
    my_abstraction_solver = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_solver"))
    start_time = time()

    # Open the file in write mode and save the data
    jldopen(filename_save, "w") do file
        file["my_abstraction_solver"] = my_abstraction_solver
    end  # This block ensures the file is closed after writing

    end_time = time()
    save_time = end_time - start_time
    @info("Time elapsed to save : $save_time")
else
    file = jldopen(filename_save, "r")
    reloaded_solver = file["my_abstraction_solver"]
    MOI.set(optimizer, MOI.RawOptimizerAttribute("abstraction_solver"), reloaded_solver)
end